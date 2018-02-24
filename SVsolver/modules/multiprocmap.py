#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parallel map with normal functional programming (without pickling functions).
# @author: "Maciek Sykulski"<macieksk@gmail.com>, around 2015

import multiprocessing


def fun(f, q_in, q_out):
    while True:
        # print("Waiting in worker...")
        i, x = q_in.get()
        if i is None:
            # print("Exiting in worker")
            q_out.put(None)
            break
        try:
            res = f(x)
        except Exception as e:
            print("Exception in Worker:{}".format(e))
        # print("Output in worker.")
        q_out.put((i, res))


def parmap(f, X, nprocs=multiprocessing.cpu_count(),
           percentage_fun=lambda i, p: None):
    q_in = multiprocessing.Queue()
    q_out = multiprocessing.Queue()

    # print("Manager starting nproc:{} processes.".format(nprocs))
    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    res = [None for i in range(len(sent))]
    for i in range(len(sent)):
        # print("Waiting in manager...(nproc:{})".format(nprocs))
        res[i] = q_out.get()
        percentage_fun(i, 100. * float(i) / len(sent))
    [q_in.put((None, None)) for _ in range(nprocs)]
    [q_out.get() for _ in range(nprocs)]
    [p.join() for p in proc]

    return [x for i, x in sorted(res)]


def test_multimap():
    def log_p(i, p):
        if i % 10 == 0:
            print("Finished {} percent".format(round(p, 2)))

    print(parmap(lambda i: i * 2, list(range(20)), percentage_fun=log_p))


def _test_treeparmap(tree, n):
    tree_nodes = list(tree.get_descendants())

    def log_p(i, p):
        if i % 10 == 0:
            print("Finished {} percent".format(round(p, 2)))

    def f(i):
        node = tree_nodes[i]
        node.add_feature("test_feature", 1)
        node.up = None
        for c in node.children:
            c.children = []
        # raise Exception("test")
        return node

    return parmap(f, list(range(n)), percentage_fun=log_p)
