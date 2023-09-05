# python3

class Edge:
    def __init__(self, u, v, capacity):
        self.u = u
        self.v = v
        self.capacity = capacity
        self.flow = 0


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
    vertex_count, edge_count = map(int, input().split())
    graph = FlowGraph(vertex_count)
    for _ in range(edge_count):
        u, v, capacity = map(int, input().split())
        graph.add_edge(u - 1, v - 1, capacity)
    return graph


def dfs(graph, u, min_capacity, visited, path):
    if u == graph.size() - 1:
        return min_capacity

    visited[u] = True
    for edge_id in graph.get_ids(u):
        edge = graph.get_edge(edge_id)
        if not visited[edge.v] and edge.capacity > edge.flow:
            path.append(edge_id)
            bottleneck = min(min_capacity, edge.capacity - edge.flow)
            result = dfs(graph, edge.v, bottleneck, visited, path)
            if result > 0:
                return result
            path.pop()

    return 0


def max_flow(graph, from_, to):
    flow = 0
    while True:
        visited = [False] * graph.size()
        path = []
        # Find an augmenting path using DFS
        augmenting_flow = dfs(graph, from_, float('inf'), visited, path)
        if augmenting_flow == 0:
            break

        # Update the flow along the augmenting path
        for edge_id in path:
            graph.add_flow(edge_id, augmenting_flow)

        flow += augmenting_flow

    return flow


if __name__ == '__main__':
    graph = read_data()
    print(max_flow(graph, 0, graph.size() - 1))
