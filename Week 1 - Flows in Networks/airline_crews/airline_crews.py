# python3

import queue


class Edge:

    def __init__(self, u, v, capacity):
        self.u = u
        self.v = v
        self.capacity = capacity
        self.flow = 0

    def available(self):
        return self.capacity - self.flow


class FlowGraph:

    def __init__(self, n):
        self.edges = []
        self.graph = [[] for _ in range(n)]

    def add_edge(self, from_, to, capacity):
        forward_edge = Edge(from_, to, capacity)
        backward_edge = Edge(to, from_, 0)
        self.graph[from_].append(len(self.edges))
        self.edges.append(forward_edge)
        self.graph[to].append(len(self.edges))
        self.edges.append(backward_edge)

    def size(self):
        return len(self.graph)

    def get_ids(self, from_):
        return self.graph[from_]

    def get_edge(self, id):
        return self.edges[id]

    def add_flow(self, id, flow):
        self.edges[id].flow += flow
        self.edges[id ^ 1].flow -= flow


def read_data():
    n, m = map(int, input().split())
    graph = FlowGraph(2 + n + m)

    for i in range(n):
        graph.add_edge(-1, i, 1)

        values = input().split()
        for idx in range(len(values)):
            if int(values[idx]) == 1:
                graph.add_edge(i, n + idx, 1)

    for i in range(m):
        graph.add_edge(n + i, n + m, 1)

    return graph, n, m


def find_path(graph, from_, to):
    visit = queue.Queue()
    visit.put((from_, []))

    visited = set()

    while not visit.empty():
        (u, p) = visit.get()
        if u in visited:
            continue

        visited.add(u)

        edges = graph.get_ids(u)

        for id in edges:
            edge = graph.get_edge(id)

            if edge.v in visited:
                continue

            if edge.available() > 0:
                if edge.v == to:
                    p.append(id)
                    return p

                next = list(p)
                next.append(id)
                visit.put((edge.v, next))

    return None


def find_min(graph, p):
    min = graph.get_edge(p[0]).available()

    for e in p:
        x = graph.get_edge(e).available()
        if x < min:
            min = x

    return min


def max_flow(graph, from_, to):
    flow = 0

    while True:
        p = find_path(graph, from_, to)
        if p is None:
            break

        min = find_min(graph, p)

        for e in p:
            graph.add_flow(e, min)

        flow += min

    return flow


if __name__ == '__main__':
    graph, n, m = read_data()
    max_flow(graph, -1, n + m)
    ans = ['-1'] * n

    for i in range(n):
        for id in graph.get_ids(i):
            edge = graph.get_edge(id)
            if edge.flow != 1 or edge.v < n:
                continue

            ans[i] = str(edge.v + 1 - n)
            break

    print(' '.join(ans))
