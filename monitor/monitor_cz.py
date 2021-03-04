#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import re
import sys
import time
import logging
import argparse
from subprocess import Popen, PIPE


def create_logger(log):
    global logger
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=log,
                        format="%(asctime)s [%(levelname)s] %(funcName)s [line: %(lineno)d]:\n    %(message)s\n",
                        filemode='w',
                        level=logging.DEBUG)
    # add a new Handler to print all INFO and above messages to stdout
    stdout = logging.StreamHandler(sys.stdout)
    stdout.setLevel(logging.INFO)
    logger.addHandler(stdout)


class TASK(object):
    '''
        job object contains basic infor, run status, and so on.

        Parameters
        ----------
        shell: string
            A string representing a job
        pattern: string
            A string identifier to indicate successfully finish.
        dependence: list
            dependent job
        run: string
            {qsub, local}

    '''
    def __init__(self, shell, pattern, dependence=None, run='qsub',
                 memory='2G', priority=None, cpu=None, queue=None):
        if not os.path.exists(shell):
            logger.error('shell file of job isn\'t exists, maybe workflow happen some error, please check!')
            sys.exit(1)
        self._shell = shell
        self.__regex__ = re.compile(pattern)    # 'pattern', finished flag in stderr file

        dependence = [] if dependence is None else dependence
        self._depend = dependence if isinstance(dependence, list) else [dependence, ]
        if run not in ('local', 'qsub'):
            logger.error('run type not in {local, qsub}, run again!')
            sys.exit(1)
        self._run = run

        self._status = 'waiting'    # (waiting, running, done, fail)
        self._jobid = None
        self._localjob = None

    @property
    def shell(self):
        return self._shell

    @property
    def dependence(self):
        return self._depend

    def add_dependence(self, depend):
        '''
            add multi dependent job for this job
        '''
        self._depend.append(depend)

    def set_dependence(self, depend):
        '''
            set dependent job from string to TASK object
        '''
        self._depend = depend

    def set_status(self, status):
        '''
            set this job running status
        '''
        assert status in ('waiting', 'running', 'done', 'fail', 'fail_dependence'), 'Error: illegal input'
        self._status = status

    def run(self):
        '''
            qsub job to calculate cluster.
        '''
        if self._run == 'qsub':
            cmd = 'qsub -o {0}.o -e {0}.e {0}'.format(self._shell)
            stdout, stderr = run_cmd(cmd).communicate()
            self._jobid = stdout.split('.')[0]
        elif self._run == 'local':
            cmd = 'sh {0} 1>{0}.o 2>{0}.e'.format(self._shell)
            job = run_cmd(cmd)
            self._jobid = job.pid
            self._localjob = job
        self.set_status('running')

    def check_status(self):
        """
            To check and set this job status.

            Returns
            -------
            job_status: string
                (waiting, running, done, fail)
        """
        if self._status == 'running':
            if self._run == 'qsub':
                cmd = 'qstat | grep {}'.format(self._jobid)
                stdout, stderr = run_cmd(cmd).communicate()
                if not stdout:
                    if self.check_signal_file():
                        self.set_status('done')
                    else:
                        self.set_status('fail')
            elif self._run == 'local':
                cmd = 'ps -u -q {0} | grep {0}'.format(self._jobid)
                stdout, stderr = run_cmd(cmd).communicate()
                if stdout.strip().endswith('<defunct>'):
                    self._localjob.wait()
                    if self.check_signal_file() and self._localjob.returncode == 0:
                        self.set_status('done')
                    else:
                        self.set_status('fail')
        return self._status

    def write_status_to_file(self, fp):
        '''
            write to status of unfinished jobs to log
        '''
        fp.write('{}\t{}\n'.format(self._shell, self._status))

    def is_finished(self):
        '''
            To check whether this job is finished, done or fail

            Returns
            -------
            flag: bool
                A boolean variable indicative of whether this job is finished.
        '''
        flag = False
        if self.check_status() in ('done', 'fail', 'fail_dependence'):
            flag = True
        return flag

    def is_normal_finished(self):
        """
            To check whether a job is normally finished.

            Returns
            -------
            flag: bool
                A boolean variable indicative of whether this job is normally finished.
        """
        flag = False
        if self.check_status() == 'done':
            flag = True
        return flag

    def satisfy_run_cond(self):
        '''
            To see whether job should be qsubbed.

            Returns
            -------
            flag: bool
                A boolean variable to indicate whether this job is to qsub.
        '''
        flag = False
        if (self._depend is None) or all([i.is_normal_finished() for i in self._depend]):
            flag = True
        return flag

    def dependence_fail(self):
        '''
            To check whether dependent jobs of this job is fail

            Returns
            -------
            flag: bool
                A boolean variable to indicate whether dependent jobs of this job is fail
        '''
        flag = False
        if self._depend:
            for i in self._depend:
                if i._status == 'fail':
                    flag = True
                    break
        return flag

    def check_signal_file(self):
        """
            Check whether error file contains finished string.

            Returns
            -------
            flag: bool
                A boolean indicative of whether error file contains finished string,
                return true if contains, otherwise false.

            Notes
            -----
            If a job has more than one stderr file, all of them will be check, this method
            will return true if one of them includes the finished string.
        """
        flag = False
        err_file = '{}.e'.format(self._shell)
        if not err_file:                # if no error file was found!
            pass
        else:
            with open(err_file) as f:
                if self.__regex__.search(f.read()):
                    flag = True
        return flag


class WORKFLOW(object):
    '''
        workflow object contains all task of project.

        Parameters
        ----------
        workflow: list
            A list representing workflow file
        project: string
            A string representing project unique ID or name
        pattern: string
            A string identifier to indicate successfully finish.
        run: string
            {qsub, local}, default: qsub
        max_job_num: int
            upper bound for the number of parallel jobs for user, default: 499
        sleeptime: int
            Seconds to sleep between two check, default: 300s
    '''
    def __init__(self, workflow, project, pattern,
                 run='qsub',
                 max_job_num=499,
                 sleeptime=300):
        self._mem = '2G'    # default memory amout
        self._pattern = pattern
        if run not in ('local', 'qsub'):
            logger.error('run type not in {local, qsub}, run again!')
            sys.exit(1)
        self._project = project
        self._max_job_num = max_job_num
        self._run = run
        self._sec = sleeptime
        self._user = self.__get_user()
        self._workflow = self.__workflow_parser(workflow)

    def __get_user(self):
        '''
        '''
        stdout, stderr = run_cmd('whoami').communicate()
        return stdout.strip()

    def __workflow_parser(self, wf):
        '''
            create a directed interdependent job list

            Returns
            -------
            workflow: list
                A list with all task object of this project.
        '''
        tasks = {}
        for i in wf:
            i = [m for m in i if m]
            if not i:
                continue
            if len(i) == 1:
                shell, memory = self.__parse_job(i[0])
                if shell not in tasks:
                    tasks[shell] = TASK(shell, self._pattern,
                                        dependence=None,
                                        run=self._run,
                                        memory=memory,
                                        priority=None,
                                        queue=None)
            elif len(i) == 2:
                shell1, memory1 = self.__parse_job(i[0])
                shell2, memory2 = self.__parse_job(i[1])
                if shell1 not in tasks:
                    tasks[shell1] = TASK(shell1, self._pattern,
                                         dependence=None,
                                         run=self._run,
                                         memory=memory1,
                                         priority=None,
                                         queue=None)
                else:
                    pass
                if shell2 not in tasks:
                    tasks[shell2] = TASK(shell2, self._pattern,
                                         dependence=shell1,
                                         run=self._run,
                                         memory=memory2,
                                         priority=None,
                                         queue=None)
                else:
                    tasks[shell2].add_dependence(shell1)
            else:
                logger.error('Workflow file format is illegal')
                sys.exit(1)

        for i in tasks:
            if tasks[i].dependence:
                # print('1.', tasks[i].shell, tasks[i].dependence)
                tasks[i].set_dependence([tasks[j] for j in tasks[i].dependence])
                # print('2.', tasks[i], tasks[i].dependence)
        return [tasks[i] for i in tasks]

    def __parse_job(self, s):
        """
            Parse a string to get the filename of the job (including absolute path) and memory.

            Parameters
            ----------
            s: string
                A string contains the job filename and associated amount memory.

            Returns
            -------
            t: tuple
                A tuple with filename (including absolute path) and memory.
        """
        memory = self._mem     # set default memory
        abs_path = s
        if s[-1] in ('G', 'g', 'M', 'm'):
            i = s.rindex(':')
            memory = s[i+1:]
            abs_path = s[:i]
        return (os.path.abspath(abs_path), memory)

    def get_number_of_running_jobs(self):
        """
            Get the total number of jobs for the user.

            Returns
            -------
            count: int
                The total number of jobs for the user.
        """
        count = 0
        if self._run == 'qsub':
            cmd = 'qstat -u {}'.format(self._user)
            stdout, stderr = run_cmd(cmd).communicate()
            return len([i for i in stdout.split('\n') if self._user in i])
        elif self._run == 'local':
            stdout, stderr = run_cmd('ps -u').communicate()
            return len(stdout.split('\n'))-10
        return count

    def is_pipeline_finished(self):
        """
            Check whether this pipeline is successfully finished.

            Returns
            -------
            flag: bool
                A boolean variable indicative of whether a pipeline is successfully finished.
        """
        flag = False
        if all([i.is_finished() for i in self._workflow]):
            flag = True
        return flag

    def monitor(self):
        """
            The main method to monitor the pipeline

            Notes
            -----
            this pipeline is normally finished or not.
        """
        # flog = open(self._project + '.task_monitor.log', 'w')
        while True:
            if self.get_number_of_running_jobs() > self._max_job_num:
                time.sleep(self._sec)
                continue

            for job in self._workflow:
                if job.satisfy_run_cond() and job.check_status() == 'waiting':
                    while True:
                        job_counts = self.get_number_of_running_jobs()
                        if job_counts >= self._max_job_num:
                            info = ("Current job counts={} exceeds maximum number allowed '{}', "
                                    "'{}' is waiting for run...\n").format(job_counts, self._max_job_num, job.shell)
                            logger.debug(info)
                            # sys.stdout.write(info); flog.write(info)
                            time.sleep(self._sec)
                        else:
                            info = "[{}]:[jod_id={}] was qsubbed.\n".format(time.asctime(), job.shell)
                            logger.info(info)
                            # sys.stdout.write(info); flog.write(info)
                            job.run()
                            break
                elif job.dependence_fail():
                    job.set_status('fail_dependence')
                    info = "[{}]:[jod_id={}] is stopped qsub, because of dependent job is fail.\n".format(time.asctime(), job.shell)
                    logger.warning(info)
                    # sys.stdout.write(info); flog.write(info)

            if self.is_pipeline_finished():
                fail_jobs = [i for i in self._workflow if i.check_status() != 'done']
                if fail_jobs:
                    info = "All jobs have finished, some jobs is fail...\n"
                    logger.warning(info)
                    # sys.stdout.write(info); flog.write(info)
                    with open("{}.unfinished_jobs.log".format(self._project), "w") as f:
                        for job in fail_jobs:
                            job.write_status_to_file(f)
                else:
                    info = "All jobs have successfully finished...\n"
                    logger.info(info)
                    # sys.stdout.write(info); flog.write(info)
                # self.send_email()
                break
            # flog.flush()
            time.sleep(self._sec)
        # flog.close()


def run_cmd(cmd):
    return Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE,
                 universal_newlines=True)    # return stdout, stderr as string.


def args_parse():
    parser = argparse.ArgumentParser(    # prog=__file__,
        usage='''
        python %(prog)s [option] -w workflow -p project''',
        description=''' submit pipeline task and monitor ''',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='''
    Note:
        e.g. the format of workflow file
            a.sh:1G    c.sh:2G
            b.sh:1G    c.sh:2G
            c.sh:3G    d.sh:2G

        This represents `c.sh` is relying on the results of `a.sh` and `b.sh`,
        and `f.sh` accepts results produced by `a.sh` and `b.sh`, a colon
        followed by is the amount of memory to feeded to associated shell script
        in invoking qsub.

        OR:
            a.sh:1G
            b.sh:2G
            c.sh:3G

        This means `a.sh`, `b.sh` and `c.sh` could be run in parallel, and no further
        shell scripts depend on neither of them.
        OR:
            a.sh:1G    c.sh:2G
            e.sh:1G
            b.sh:1G    c.sh:2G
            f.sh:1G
            c.sh:3G    d.sh:2G

        This indicates 'c.sh' depends upon 'a.sh' and 'b.sh' while 'd.sh' relies on
        'c.sh', and none of the script depends on 'e.sh' and 'f.sh', so 'e.sh', 'f.sh'
        can be run in parallel with 'a.sh' and 'b.sh'. ''')
    parser.add_argument('-w', '--workflow',
                        dest='workflow', metavar='',
                        help='analysis task workflow file, [REQUIRED]')
    parser.add_argument('-p', '--project',
                        dest='project', metavar='',
                        help='project name for ouput log, [REQUIRED]')
    parser.add_argument('--pattern',
                        dest='pattern', metavar='',
                        default='Across the Great Wall we can reach every corner in the world',
                        help='job finished string, Default: Across the Great Wall we can reach every corner in the world')
    parser.add_argument('--run',
                        dest='run', metavar='',
                        choices=('qsub', 'local'), default='local',
                        help='set running type, select from {qsub, local}, default: local')
    parser.add_argument('--sleep',
                        dest='sleep', metavar='',
                        type=int, default=300,
                        help='sleeping interval in sec, default: 300')
    parser.add_argument('--max_job',
                        dest='maxjob', metavar='',
                        type=int, default=499,
                        help='maximum number of jobs allowed, default: 499')
    args = parser.parse_args()
    if not all([args.workflow, args.project]):
        parser.print_help()
        sys.exit(1)
    return args.workflow, args.project, args.pattern, args.run, args.sleep, args.maxjob


def main():
    workflow, project, pattern, run, sleep, maxjob = args_parse()
    create_logger('{}.log'.format(project))

    if not os.path.exists(workflow):
        logger.error('Workflow file isn\'t exists, please check!')
        sys.exit(1)
    with open(workflow) as f:
        workflow = [i.strip().split() for i in f if (i.strip() and not i.startswith('#'))]

    WF = WORKFLOW(workflow, project, pattern,
                  run=run,
                  max_job_num=maxjob,
                  sleeptime=sleep)
    WF.monitor()


if __name__ == '__main__':
    main()
