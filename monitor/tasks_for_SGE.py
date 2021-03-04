#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
from os.path import basename, splitext, exists, abspath


class tasks(object):
    '''
    '''
    def __ini__(self, cmd):
        self.cmd = cmd

    def cmd_execute(self, wait=True):
        try:
            assert not os.system(self.cmd)
        except Exception as e:
            info = ('cmd execute fail... \n'
                    '    cmd: "{}" \n'
                    '    err: "{}" \n').format(self.cmd, e)
            raise Exception(info)


def valid_workflows(workflows):
    '''
    workflows:    list
        workflow list, [workflow, workflow, ...]

        workflow:
            [shell1, [shell2, shell3, ...], ...]

    return workflows
    '''
    if not isinstance(workflows, list):
        raise Exception('Unsupport workflows type, must be list in workflows')

    def valid_workflow(order_tasks):
        if not set([type(i) for i in order_tasks]).issubset({list, str}):
            raise Exception(('Unsupport workflow type, '
                            'must be list or str in workflow')
                            )
        tmp = [i for i in order_tasks if isinstance(i, list)]
        if tmp and not ({str} == set([type(j) for i in tmp for j in i])):
            raise Exception('Unsupport tasks type, must be str in tasks')
        return True

    if {str}.issubset(set([type(i) for i in workflows])):
        if valid_workflow(workflows):
            return [workflows, ]
    else:
        if all([valid_workflow(workflow) for workflow in workflows]):
            return workflows


def save_order_task(workflows, output, order_main=['A', ]):
    '''
    workflows:    list
        workflow list, [workflow, workflow, ...]

        workflow:
            [shell1, [shell2, shell3, ...], ...]

    output:    string
        output file of order task
        shell | order | order_main
    '''
    workflows = valid_workflows(workflows)
    assert len(workflows) == len(order_main), \
            'Inconsistent length of "workflows" and "order_main"'

    def single_workflow(tasks, main):
        ordertasks = []
        for index, task in enumerate(tasks):
            if isinstance(task, str):
                ordertasks.append('{} | {} | {}'.format(task, index, main))
            elif isinstance(task, list):
                for shell in task:
                    ordertasks.append('{} | {} | {}'.format(shell, index, main))
        return ordertasks

    workflow_txt = []
    for index, main in enumerate(order_main):
        workflow = workflows[index]
        workflow_txt += single_workflow(workflow, main)
    with open(output, 'w') as f:
        f.write('\n'.join(workflow_txt)+'\n')


def read_order_task(workflowsfile):
    '''
    workflowsfile:    string
        file of order task, from Function <save_order_task> output
        file format:
            shell1 | order1 | order_main1
            shell2 | order2 | order_main2

    return:    list
        workflows list, [workflow, workflow, ...]

        workflow:
            [shell1, [shell2, shell3, ...], ...]
    '''
    from collections import defaultdict

    def single_workflow(order_task):
        '''
        order_task:    list
            task list, [[shell1,order1,order_main1],
                        [shell2,order2,order_main1], ...]
        return
        '''
        workflow = defaultdict(list)
        for task in order_task:
            shell, order = task[:2]
            workflow[order].append(shell)
        workflow = [i[1] for i in sorted(workflow.items(),
                                         key=lambda x: int(x[1]))]
        for index, task in enumerate(workflow):
            if len(task) == 1:
                workflow[index] = task[0]
        return workflow

    with open(workflowsfile) as f:
        tasks = [task.split('|') for task in f if task.strip() != '']
        tasks = map(lambda x: [i.strip() for i in x], tasks)

        workflow_multi = defaultdict(list)
        for task in tasks:
            workflow_multi[task[2]].append(task)

        return [single_workflow(workflow_multi[i]) for i in workflow_multi]


class SGE_task(object):
    """
    name:    string
        name of task, must be unique in all tasks name
    cmd:    string
        command line must be start with parser program, sh, python, perl et al.
    """
    def __init__(self, name, cmd):
        self.name = name
        if cmd.lower().endswith('.sh'):
            self.cmd = 'sh {}'.format(cmd)
        elif cmd.lower().endswith('.py'):
            self.cmd = 'python {}'.format(cmd)
        elif cmd.lower().endswith('.pl'):
            self.cmd = 'perl {}'.format(cmd)
        elif cmd.lower().endswith('.r'):
            self.cmd = 'Rscript {}'.format(cmd)
        else:
            raise Exception('Can\'t find parser for CMD:\n  {}!\n'.format(cmd))

    def sge(self):
        txt = ('job_begin\n'
               '    name {}\n'
               '    sched_options -cwd -V -q all.q \n'
               '    cmd {} \n'
               'job_end\n').format(self.name, self.cmd)
        return txt


def OrderWorkflow4SGE(workflow):
    '''
    convert tasks list to sge order format

    workflow:    list
        [shell1, [shell2, shell3, ...], ...]

    return:    tupe
        (sge_task, sge_order)
        sge_task:   list
        sge_order:  list
    '''
    workflow = valid_workflows([workflow, ])[0]
    TaskName = lambda x: basename(splitext(x)[0])
    # convert task list of workflow to SGE task list
    sge_task = []
    for index, task in enumerate(workflow):
        if isinstance(task, str):
            tmp = SGE_task(TaskName(task), task)
            workflow[index] = tmp
            sge_task.append(tmp.sge())
        elif isinstance(task, list):
            tmp = [SGE_task(TaskName(i), i) for i in task]
            workflow[index] = tmp
            sge_task += [i.sge() for i in tmp]
        else:
            raise Exception('Unsupport input type of Task ')

    # order SGE task
    SmapleName = lambda name: name.partition('_')[2]
    SGEOrder = lambda first,later: 'order {} after {}\n'.format(later, first)

    def links_taskflows(upstream, downstream):
        fake_task = ('job_begin\n'
                     '    name Fake_task\n'
                     '    host localhost\n'
                     '    cmd echo Shit is comming...\n'
                     'job_end\n')
        fake_name = 'Fake_task'
        OrderList = [SGEOrder(i.name, fake_name) for i in upstream]
        OrderList += [SGEOrder(fake_name, i.name) for i in downstream]
        return OrderList, fake_task

    sge_order = []
    startup = workflow[0]
    for task in workflow[1:]:
        if isinstance(startup, list) and isinstance(task, str):
            sge_order += [SGEOrder(i.name, task.name) for i in startup]
        elif isinstance(startup, str) and isinstance(task, list):
            sge_order += [SGEOrder(startup.name, i.name) for i in task]
        elif isinstance(startup, str) and isinstance(task, str):
            sge_order.append(SGEOrder(startup.name, task.name))
        elif isinstance(startup, list) and isinstance(task, list):
            if len(startup) == len(task):
                startup_dict = {SmapleName(i.name): i.name for i in startup}
                task_dict = {SmapleName(i.name): i.name for i in task}
                if set(startup_dict.keys()) == set(task_dict.keys()):
                    sge_order += [SGEOrder(startup_dict[i], task_dict[i])
                                  for i in startup_dict]
                else:
                    OrderList, fake_task = links_taskflows(startup, task)
                    sge_order += OrderList
                    sge_task.append(fake_task)
            else:
                OrderList, fake_task = links_taskflows(startup, task)
                sge_order += OrderList
                sge_task.append(fake_task)

        startup = task
    return sge_task, sge_order


def OutSGEOrder2File(workflows, jobs, log, shell=None):
    '''
    write "workflows" to "jobs", which is available for SGE Job Manager

    workflows:    list
        [workflow, workflow, ...]

        workflow:    list
            [shell1, [shell2, shell3, ...], ...]

    jobs:    string
        output file,

    log:    string
        output logging dir for analysis Process
    '''
    if not exists(log):
        os.mkdir(log)

    workflows = valid_workflows(workflows)
    TasksDetail, TasksOrder = [], []
    for workflow in workflows:
        task, task_order = OrderWorkflow4SGE(workflow)
        TasksDetail += task
        TasksOrder += task_order

    tmp = []
    for i in TasksDetail:
        if i not in tmp:
            tmp.append(i)
    TasksDetail = '\n'.join(tmp)

    tmp = []
    for i in TasksOrder:
        if i not in tmp:
            tmp.append(i)
    TasksOrder = ''.join(tmp)

    log = 'log_dir {}'.format(log)

    with open(jobs, 'w') as f:
        f.write(('{task_detail}\n\n'
                 '{task_order}\n'
                 '{logdir}\n').format(task_detail=TasksDetail,
                                      task_order=TasksOrder,
                                      logdir=log)
                )
    if shell:
        with open(shell, 'w') as f:
            f.write('sjm {jobs}'.format(**locals()))
        assert not os.system('chmod +x {shell}'.format(**locals()))


def main():
    import argparse
    parser = argparse.ArgumentParser(
                                    formatter_class=argparse.RawTextHelpFormatter,
                                    epilog='''
    #################################################################
    workflows file:
        shell1 | order1 | order_mainA
        shell2 | order2 | order_mainA
        shell3 | order3 | order_mainA
        shell4 | order4 | order_mainB

      shell:      string, linux shell file
      order:      int, sort task in same analysis block
      order_main: string or int, distinguish different analysis block
    ################################################################# '''
                                    )
    parser.add_argument('-p', '--project',
                        dest='project',
                        metavar='',
                        required=True,
                        help='Unique project ID, [REQUIRED]'
                        )
    parser.add_argument('-w', '--workflows',
                        dest='workflows',
                        metavar='',
                        required=True,
                        help='workflows file, [REQUIRED]'
                        )
    parser.add_argument('--root',
                        dest='root_dir',
                        metavar='',
                        default='./',
                        help='root dir of analysis, default: current path'
                        )
    parser.add_argument('--log',
                        dest='log_dir',
                        metavar='',
                        help='the dir for output analysis logging'
                        )
    select = parser.add_argument_group('Select Task run type')
    RunSelect = select.add_mutually_exclusive_group()
    RunSelect.add_argument('--sge',
                           dest='sge',
                           action='store_true',
                           help='use SGE submit tasks, [Default]'
                           )
    RunSelect.add_argument('--locate',
                           dest='locate',
                           action='store_true',
                           help='use locate type run task, Function uncomplete'
                           )
    args = parser.parse_args()
    project = args.project
    workflowsfile = args.workflows
    root_dir = abspath(args.root_dir)
    if not exists(root_dir):
        os.mkdir(root_dir)
    log_dir = args.log_dir if args.log_dir else '{}/log'.format(root_dir)
    if not exists(log_dir):
        os.mkdir(log_dir)
    sge = args.sge
    locate = args.locate
    if not all([sge, locate]):
        sge = True

    if sge:
        workflows = read_order_task(workflowsfile)
        OutSGEOrder2File(workflows,
                         '{}/{}_analysis.JOB'.format(log_dir, project),
                         log_dir,
                         shell='{}/{}_analysis.sh'.format(root_dir, project)
                         )
    elif locate:
        sys.stderr.write('Function uncomplete, Waiting for add...')
        sys.exit(1)

if __name__ == '__main__':
    main()
