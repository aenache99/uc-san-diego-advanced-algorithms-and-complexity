# python3
class StockCharts:
    def read_data(self):
        n, k = map(int, input().split())
        stock_data = [list(map(int, input().split())) for i in range(n)]
        return stock_data

    def write_response(self, result):
        print(result)

    def can_place_above(self, stock1, stock2):
        # Check if stock1 can be placed above stock2 without intersecting
        return all(x > y for x, y in zip(stock1, stock2))

    def build_bipartite_graph(self, stock_data):
        n = len(stock_data)
        m = len(stock_data[0])
        graph = [[False] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                if self.can_place_above(stock_data[i], stock_data[j]):
                    graph[i][j] = True

        return graph

    def max_bipartite_matching(self, graph):
        n = len(graph)
        m = len(graph[0])
        matching = [-1] * m
        visited = [False] * n

        def dfs(u):
            visited[u] = True
            for v in range(m):
                if graph[u][v] and (matching[v] == -1 or (not visited[matching[v]] and dfs(matching[v]))):
                    matching[v] = u
                    return True
            return False

        for u in range(n):
            visited = [False] * n
            dfs(u)

        return matching.count(-1)

    def min_charts(self, stock_data):
        graph = self.build_bipartite_graph(stock_data)
        max_matching = self.max_bipartite_matching(graph)
        return max_matching

    def solve(self):
        stock_data = self.read_data()
        result = self.min_charts(stock_data)
        self.write_response(result)


if __name__ == '__main__':
    stock_charts = StockCharts()
    stock_charts.solve()
