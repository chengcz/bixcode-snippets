#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function, division
import os
from os.path import dirname, basename, exists, abspath, isfile
import sys
import json
import yaml    #
import time
import hashlib
import atexit
import signal
import logging
import subprocess
from graphviz import Digraph, Graph

try:
    import cPickle as pickle
except:
    import pickle

# import transaction
# from persistent import Persistent
# from ZODB import FileStorage, DB
# from apscheduler.schedulers.background import BackgroundScheduler
#



def create_logger(log):
    global logger
    fmt='%(levelname)-5s @ %(asctime)s, %(filename)s:%(lineno)s, %(message)s',
    logging.basicConfig(
        filename=log,
        filemode='a',
        format=fmt,
        level=logging.DEBUG,
        # stream=sys.stderr
    )
    logger = logging.getLogger(__name__)
    # streamlog = logging.StreamHandler(sys.stderr)
    # # streamlog = logging.FileHandler(log)
    # streamlog.setLevel(logging.INFO)
    # streamlog.setFormatter(logging.Formatter(fmt))
    # logger.addHandler(streamlog)


fmt='%(levelname)-5s @ %(asctime)s, %(filename)s:%(lineno)s, %(message)s',
logging.basicConfig(
    format=fmt,
    level=logging.DEBUG,
    stream=sys.stderr
)
logger = logging.getLogger(__name__)


class pickleDB(obejct):
    def __init__(self):
        pass

    def write(self, obj, files):
        with open(files, 'wb') as fo:
            return pickle.dump(obj, fo)

    def read(self, files):
        with open(files, 'rb') as fi:
            return pickle.load(fi)


class ZODatabase(object):
	'''
    ref
    ---
    https://www.ibm.com/developerworks/cn/aix/library/au-zodb/index.html
    '''
	def Init(self, filename):
		try:
			self.storage = FileStorage.FileStorage(filename)
		except:
			return 0
		else:
			self.db = DB(self.storage)
			return 1

	def Open(self):
		self.conn = self.db.open()
		return self.conn.root()

	def Close(self):
		transaction.commit()
		self.conn.close()

	def Clean(self):
		self.db.pack()

	def Uninit(self):
		self.storage.close()


class Task(object):
    '''
    parameter
    ---------
    shell:     [string], command line / shell file, [Required]
    input:     [list/tuple/string], input files of CMD
    output:    [list/tuple/string], output files of CMD
    '''

    slots = ('_shell', '_input', '_output', '_depend','_attr')
    def __init__(self, shell, **kwargs):
        self._shell = shell
        self._input = kwargs.pop('input', ())
        self._output = kwargs.pop('output', ())
        self._depend = ()
        self._attr = kwargs
        self.__parameter_validation()

    def __parameter_validation(self):
        if not isinstance(self._shell, str):
            raise TypeError('Only support string for Task.shell ')
        if isfile(self._shell):
            self._attr['shell'] = True
            self._attr['name'] = self._attr.get('name', self._shell)
            with open(self._shell) as fi:
                cmd = fi.read().strip()
            self._attr['id'] = self.__md5sum(cmd)
        else:
            self._attr['shell'] = False
            self._attr['id'] = self.__md5sum(self._shell)
        self._input = self.__dtype_convertion(self._input)
        self._output = self.__dtype_convertion(self._output)

    @staticmethod
    def __dtype_convertion(idat):
        if not isinstance(idat, (str, tuple, list)):
            msg = 'Only support list/tuple/string date for cmd Task.input / Task.output.'
            raise TypeError(msg)
        if isinstance(idat, str):
            idat = idat.split(',')
        return tuple(set(idat))

    @staticmethod
    def __md5sum(string):
        if isinstance(string, str): string = string.encode('utf-8')
        md5value = hashlib.md5(string).hexdigest()
        return md5value

    def __repr__(self):
        msg = f"<monitor.Task object, Name: {self.name}>"
        return msg

    def __str__(self):
        msg = (
            "<monitor.Task object>\n"
            f"    task ID: {self.id}\n"
            f"    Name:    {self.name}\n"
            f"    cmd:     {self.shell}\n"
            f"    status:  {self.status}\n"
        )
        return msg

    @property
    def id(self):
        return self._attr['id']

    @property
    def name(self):
        return self._attr.get('name', self.id)

    @property
    def shell(self):
        return self._shell

    @property
    def input(self):
        return self._input

    @property
    def output(self):
        return self._output

    @property
    def depend(self):
        return self._depend

    @property
    def status(self):
        return self._attr.get('status', 'init')

    @property
    def rid(self):
        return self._attr.get('runid', None)
    
    @property
    def time(self):
        return self._attr.get('time', (None, None))
    
    @property
    def priority(self):
        return self._attr.get('priority', 0)

    @staticmethod
    def __append_cmd_files(source, files):
        if isinstance(files, (str, list, tuple, set)):
            if isinstance(files, str):
                files = files.split(',')
            return tuple(set(source) | set(files))
        else:
            raise TypeError('Not support data type for Task.input / Task.ouput')

    def add_input(self, files):
        self._input = self.__append_cmd_files(self._input, files)

    def add_output(self, files):
        self._output = self.__append_cmd_files(self._output, files)

    def add_depend(self, jobs):
        '''
        jobs:    single Task object or Task list
        '''
        if isinstance(jobs, Task):
            jobs = [jobs, ]
        elif isinstance(jobs, (list, set, tuple)):
            for job in jobs:
                assert isinstance(job, Task)
        else:
            raise TypeError('Not support data type for Task.depend')
        tmp = set(self._depend)
        tmp.update(jobs)
        self._depend = tuple(tmp)

    def to_json(self):
        task = self._attr
        task['input'] = self.input
        task['output'] = self.output
        task['shell'] = self.shell
        task['depend'] = [x.shell for x in self.depend]
        task['name'] = self.name
        return task

    def run(self):
        pass

    def submited(self):
        pass

    def satisfy_run_cond(self):
        pass

    def is_finished(self):
        pass

    def is_normal_finished(self):
        pass

    def set_status(self):
        pass

    def check_status(self):
        pass


class Workflow(object):
    '''

    parameter
    ---------
    project:   string
    tasks:     tasks file, yaml, json, ordered shell
    workdir:   workdir for outpout
    priority:  intger

    ref
    ---
    network:    https://www.jianshu.com/p/e543dc63454f
    '''
    __slots__ = ('_project', '_tasks', '_workdir', '_attr')
    def __init__(self, project, tasks, workdir, **kwargs):
        self._project = project
        self._tasks = self.__parse_tasks_from_tasksfile(tasks)
        self._workdir = abspath(workdir)
        self._attr = kwargs
        self.__space_prepare()

    def __parse_tasks_from_tasksfile(self, tasksfile):
        '''
        input style of tasks file
        -------------------------
        a. BGI txt style, "task shell1    task shell2"
        b. yaml/json style, Use input/ouput file to define task order
        '''
        if tasksfile.lower().endswith('txt'):
            return self.__Ordered_Tasks(tasksfile)
        elif tasksfile.lower().endswith(('yaml', 'json')):
            return self.__Input_Output_Tasks(tasksfile)
        else:
            raise TypeError('Not support file format(suffix) for Workflow.tasks')

    def __space_prepare(self):
        if not isinstance(self._project, str):
            raise TypeError('')
        self.__makedirs(self._workdir)
        self._attr['logdir'] = f'{self._workdir}/log'
        self.__makedirs(self._attr['logdir'])
        self._attr['log'] = f'{self._workdir}/log/{self._project}.log'

    @staticmethod
    def __makedirs(path):
        try:
            os.makedirs(path, mode=0o755, exist_ok=True)
        except Exception as e:
            logger.error('{}, now exit!'.format(e))
            sys.exit(1)   

    def __repr__(self):
        msg = f"<monitor.Workflow object, project: {self.project}>"
        return msg

    def __str__(self):
        msg = (
            "<monitor.Workflow object>\n"
            f"    project:      {self._project}\n"
            f"    task number:  {len(self._tasks)}\n"
            f"    workdir:      {self._workdir}\n"
        )
        return msg

    def __iter__(self):
        for job in self._tasks:
            yield job

    @property
    def project(self):
        return self._project

    @property
    def workdir(self):
        return self._workdir

    @property
    def tasks(self):
        return self._tasks
    
    @property
    def logdir(self):
        return self._attr['logdir']

    @property
    def priority(self):
        return self._attr.get('priority', 0)

    @staticmethod
    def __Ordered_Tasks(tasksfile):
        with open(tasksfile) as fp:
            ordered = [x.strip().split('\t') for x in fp if not x.startswith('#')]
            tasks = set([y for x in ordered for y in x])
        tkdepent = {}
        for depentask in ordered:
            fir, sec = map(lambda x: x.partition(':')[0], depentask[:2])
            tkdepent.setdefault(sec, []).append(fir)
        tkCollapse = {}
        for task in tasks:
            shell, _, memory = task.partition(':')
            tkCollapse[shell] = Task(shell, name=basename(shell), memory=memory)
        iniTasks = []
        for name, task in tkCollapse.items():
            depent = [tkCollapse[x] for x in tkdepent.get(name, ())]
            task.add_depend(depent)
            iniTasks.append(task)
        return tuple(iniTasks)

    @staticmethod
    def __Input_Output_Tasks(tasksfile):
        '''
        Use input/ouput file to define task order
        taks1: {input: xxx, output: xxx, shell: xxx, ... }
        '''
        with open(tasksfile) as fp:
            if tasksfile.lower().endswith('yaml'):
                tasks = yaml.load(fp)
            elif tasksfile.lower().endswith('json'):
                tasks = json.load(fp)
        iniTasks = [Task(tasks[x].pop('shell'), **tasks[x]) for x in tasks]
        tkoutput = {}, {}
        for task in iniTasks:
            for m in task.output:
                tkoutput.setdefault(m, []).append(task)
        for task in iniTasks:
            depend = [y for x in task.input for y in tkoutput[x]]
            task.add_depend(depend)
        return tuple(iniTasks)

    def flow_chart(self):
        '''
        ref
        ---
        https://www.ibm.com/developerworks/cn/aix/library/au-aix-graphviz/index.html
        '''
        dot = Digraph(comment=self.project)
        colors = {
            'init':    '',
            'pending': '#ff9800',
            'running': '#03a9f4',
            'fail':    '#f44336',
            'success': '#449d48',
        }
        for index, job in enumerate(self._tasks):
            # job,    Task class
            dot.node(
                job.id, job.name, # shape='', color='',
                style='filled', fillcolor=colors.get(job.status, '#90908b')
            )
            for dpjob in job.depend:
                dot.edge(dpjob.id, job.id, constraint='true')
        sys.stdout.write(dot.source)
        # flowfile = f'{self.workdir}/{self.project}.dot'
        # cmd = (f'python %(prog)s monitor -p {self.project} --dag | '
        #         'dot –Tsvg –o workflow.svg')
        # subprocess.check_call(cmd, shell=True)

    def stats(self, output=False):
        sts = {}
        for job in self._tasks:
            sts[job.status] = sts.get(job.status, 0) + 1
        if output:
            sys.stdout.write(f'{"Status": <15}\t{"Count": <10}\n')
            for status, count in sorted(sts.items(), key=lambda x: x[0]):
                sys.stdout.write(f'{status: <15}\t{count: <10}\n')
        return sts

    def to_json(self):
        tasks = {x.name: x.to_json() for x in self._tasks}
        return tasks

    def logger(self):
        pass

    def kill(self):
        pass

    def start(self):
        pass

    def stop(self):
        pass

    def pause(self):
        pass


class Submiter(object):
    '''

    local:
    SGE:     https://drmaa-python.readthedocs.io/en/latest/tutorials.html
    Slurm:
    '''
    pass


class Daemon(object):
    def __init__(self, pidlock,
        stdin='/dev/null', stdout='/dev/null', stderr='/dev/null'
    ):
        self._lockfile = pidlock
        self._stdin = stdin
        self._stdout = stdout
        self._stderr = stderr

    def Daemonize(self):
        if exists(self._lockfile):
            raise RuntimeError('This project daemon is running.')
        # First fork (detaches from parent)
        try:
            if os.fork() > 0:
                raise SystemExit(0)
        except OSError as e:
            raise RuntimeError('fork #1 faild: {0} ({1})'.format(e.errno, e.strerror))
        os.chdir('/')
        os.setsid()
        os.umask(0o22)
        # Second fork (relinquish session leadership)
        try:
            if os.fork() > 0:
                raise SystemExit(0)
        except OSError as e:
            raise RuntimeError('fork #2 faild: {0} ({1})\n'.format(e.errno, e.strerror))
        # Flush I/O buffers
        sys.stdout.flush()
        sys.stderr.flush()
        # Replace file descriptors for stdin, stdout, and stderr
        with open(self._stdin, 'rb', 0) as fp:
            os.dup2(fp.fileno(), sys.stdin.fileno())
        with open(self._stdout, 'ab', 0) as fp:
            os.dup2(fp.fileno(), sys.stdout.fileno())
        with open(self._stderr, 'ab', 0) as fp:
            os.dup2(fp.fileno(), sys.stderr.fileno())
        # Write the PID file
        with open(self._lockfile, 'w') as fo:
            fo.write(os.getpid())
        # Arrange to have the PID file removed on exit/signal
        atexit.register(lambda: os.remove(self._lockfile))
        signal.signal(signal.SIGTERM, self.__sigterm_handler)

    # Signal handler for termination (required)
    @staticmethod
    def __sigterm_handler(signo, frame):
        raise SystemExit(1)


class scheduler(Daemon):
    def __init__(self):
        pass

    


class monitor(Daemon):
    def __init__(self):
        self._sched = self.__Scheduler()
        # daemonize
        self._lockfile = '/tmp/monitor.pid'.format(project)
        self._stdin = stdin
        self._stdout = stdout
        self._stderr = stderr

    def __Scheduler(self):
        sched = BackgroundScheduler()
        sched.add_job(self.monitor, 'interval', minutes=30, id=self.project)
        return sched

    def logger(self):
        pass

    def number_of_running_jobs(self):
        pass

    def run(self):
        self._sched.start()

    def start(self):
        try:
            self.Daemonize()
        except RuntimeError as e:
            logger.error(e)
            raise SystemExit(1)
        self.run()

    def stop(self, waiting=True):
        '''
        '''
        self._sched.remove_job(self.project)
        if exists(self._lockfile):
            with open(self._lockfile) as f:
                pid = int(f.read())
            try:
                os.kill(pid, signal.SIGTERM)
            except OSError as e:
                logger.warning(e)
                os.remove(self._lockfile)
        else:
            logger.info('Not running, exit now.')
            raise SystemExit(1)

    def pause(self):
        self._sched.pause_job(self.project)

    def resume(self):
        self._sched.resume_job(self.project)

    def restart(self):
        self.stop()
        self.start()

    def email(self):
        pass


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument()
