# python3
from time import sleep

DEBUG = False


class AdAllocation:
    def __init__(self, A, b, c):
        self.A = A
        self.b = b
        self.c = c

    def solve(self):
        n = len(self.A)
        m = len(self.A[0])
        sign = ["<=" for _ in range(n)]
        nonneg_x = [True for _ in range(m)]
        maximize = True
        algo = SIMPLEX(n, m, self.A, sign, self.b, nonneg_x, self.c, maximize)
        res = algo.solve()
        if res[1] == "No solution":
            print("No solution")
        elif res[1] == "Infinity":
            print("Infinity")
        else:
            print("Bounded solution")
            print(" ".join(map(str, res[0][1])))


class SIMPLEX:
    EPS = 0.00001

    def __init__(self, n, m, A, sign, b, nonneg_x, c, maximize=True):
        # the number of inequalities
        self.n = n
        # the number of variables
        self.m = m
        # coefficients of linear inequalities
        # Ax <= b, where A (nxm), x (m), b (n)
        self.A = A
        self.sign = sign
        self.b = b
        # non-negativity constraints on variables
        self.nonneg_x = nonneg_x
        # the coefficients of the objective to maximize
        # c (m)
        self.c = c
        self.maximize = maximize

    def solve(self):
        res = None

        A, b, c = self._to_standard_form(
            self.A,
            self.sign,
            self.b,
            self.nonneg_x,
            self.c,
            self.maximize
        )
        self.print_system_standard(A, b, c)

        A, b, c, N, B, v, msg = self._init_simplex(A, b, c)

        if msg == "No solution":
            return res, msg

        while True:
            if DEBUG: self.print_canonical_form(A, b, c, N, B, v)

            entering_var = self._get_entering_var(c, N)
            if entering_var is None:
                # optimization finished
                x = [0 for _ in range(len(A) + len(A[0]))]
                for i, j in enumerate(B):
                    x[j] = b[i]

                res = (v, x[:len(A[0])])
                msg = "Bounded solution"

                break

            if DEBUG: print(f"Entering: {entering_var + 1}")

            leaving_var = self._get_leaving_var(A, b, N, B, entering_var)
            if leaving_var is None:
                # optimization is unconstrained
                res = float("inf")
                msg = "Infinity"
                break
            if DEBUG: print(f"Leaving: {leaving_var + 1}")

            A, b, c, v, N, B = self._pivot(A, b, c, v, N, B, leaving_var, entering_var)

            if DEBUG: sleep(0.1)

        return res, msg

    @staticmethod
    def print_system_standard(A, b, c):
        if DEBUG: print()
        obj_msg = "maximize "
        for i, c_i in enumerate(c):
            if i > 0:
                if c_i >= 0:
                    obj_msg += " + "
                else:
                    obj_msg += " - "
                    c_i *= -1
            obj_msg += f"{c_i}*x{i + 1}"
        if DEBUG: print(obj_msg)

        for i, row in enumerate(A):
            row_msg = ""
            for j, x_i in enumerate(row):
                if j > 0:
                    if x_i >= 0:
                        row_msg += " + "
                    else:
                        row_msg += " - "
                        x_i *= -1
                row_msg += f"{x_i}*x{j + 1}"
            row_msg += " <= "
            row_msg += f"{b[i]}"
            if DEBUG: print(row_msg)

    @staticmethod
    def print_canonical_form(A, b, c, N, B, v):
        if DEBUG: print()
        obj_msg = f"z = {v} "
        for i, c_i in enumerate(c):
            if c_i >= 0:
                obj_msg += " + "
            else:
                obj_msg += " - "
                c_i *= -1
            obj_msg += f"{c_i}*x{N[i] + 1}"
        if DEBUG: print(obj_msg)

        for i, row in enumerate(A):
            row_msg = f"x{B[i] + 1} = {b[i]} - ("
            for j, x_i in enumerate(row):
                if j > 0:
                    if x_i >= 0:
                        row_msg += " + "
                    else:
                        row_msg += " - "
                        x_i *= -1
                row_msg += f"{x_i}*x{N[j] + 1}"
            row_msg += ")"
            if DEBUG: print(row_msg)

    @staticmethod
    def _to_standard_form(A, sign, b, nonneg_x, c, maximize):
        A = A.copy()
        b = b.copy()
        c = c.copy()

        # 1. Objective function must be maximized.
        if not maximize:
            c = [x * (-1) for x in c]

        # 2. All variables must have non-negativity inequalities.
        # Replace each decision variable unconstrained in sign by
        # a difference between two non-negative variables.
        # This replacement applies to all equations including the
        # objective function.
        shift = 0
        for i, nonneg in enumerate(nonneg_x):
            if not nonneg:
                j = i + shift
                shift += 1

                for k, row in enumerate(A):
                    A[k] = row[:(j + 1)] + [(-1) * row[j]] + row[(j + 1):]

                c = c[:(j + 1)] + [(-1) * c[j]] + c[(j + 1):]

        # 3. Replace equalities with inequalities.
        sign_ = []
        A_ = []
        b_ = []
        for i, s in enumerate(sign):
            if s == "=":
                sign_.append("<=")
                A_.append(A[i])
                b_.append(b[i])

                sign_.append(">=")
                A_.append([(-1) * x for x in A[i]])
                b_.append((-1) * b[i])
            else:
                sign_.append(s)
                A_.append(A[i])
                b_.append(b[i])
        sign = sign_
        A = A_
        b = b_

        # 4. All constraints must be inequalities 'less than or equal'.
        for i, s in enumerate(sign):
            if sign == ">=":
                A[i] = [(-1) * x for x in A[i]]
                b[i] *= -1

        return A, b, c

    def _init_simplex(self, A, b, c):
        n = len(A)
        m = len(A[0])

        min_b = float("inf")
        k = None
        for i, b_i in enumerate(b):
            if b_i < min_b:
                min_b = b_i
                k = i

        if min_b >= 0:
            msg = "System has solution"

            N = list(range(m))
            B = list(range(m, m + n))
            v = 0

            return A, b, c, N, B, v, msg

        if DEBUG: print("START INIT PROCEDURE")

        c_ = [0 for _ in range(m)]
        c_.append(-1)

        for row in A:
            row.append(-1)
        B = list(range(m + 1, m + 1 + n))
        N = list(range(1, m + 1))
        N.append(0)
        v = 0
        leaving_var = m + k + 1
        entering_var = 0

        if DEBUG: self.print_canonical_form(A, b, c_, N, B, v)
        A, b, c_, v, N, B = self._pivot(A, b, c_, v, N, B, leaving_var, entering_var)

        while True:
            if DEBUG: self.print_canonical_form(A, b, c_, N, B, v)

            entering_var = self._get_entering_var(c_, N)
            if entering_var is None:
                # optimization finished
                res = v
                msg = "Bounded solution"
                break
            if DEBUG: print(f"Entering: {entering_var + 1}")

            leaving_var = self._get_leaving_var(A, b, N, B, entering_var)
            if leaving_var is None:
                # optimization is unconstrained
                res = float("inf")
                msg = "Infinity"
                break
            if DEBUG: print(f"Leaving: {leaving_var + 1}")

            A, b, c_, v, N, B = self._pivot(A, b, c_, v, N, B, leaving_var, entering_var)

            if DEBUG: sleep(0.1)

        x = [0 for _ in range(n + m + 1)]
        for i, j in enumerate(B):
            x[j] = b[i]

        if -self.EPS <= x[0] <= self.EPS:
            if 0 in B:
                if DEBUG: print("Final pivot")
                entering_var = -1
                row_i = B.index(0)
                for i in range(m + 1):
                    if A[row_i][i] != 0:
                        entering_var = N[i]
                        break
                A, b, c_, v, N, B = self._pivot(A, b, c_, v, N, B, 0, entering_var)

            msg = "System has solution"
            aux_var_col_i = N.index(0)
            N.remove(0)
            N = [x - 1 for x in N]
            B = [x - 1 for x in B]
            for i in range(len(A)):
                A[i] = A[i][:aux_var_col_i] + A[i][aux_var_col_i + 1:]
            c_ = [0 for _ in range(m)]
            v = 0
            for i, c_i in enumerate(c):
                if c_i != 0:
                    if i in N:
                        c_[N.index(i)] += c_i
                    else:
                        v += c_i * b[B.index(i)]
                        for j in range(m):
                            c_[j] += -1 * c_i * A[B.index(i)][j]

            if DEBUG: print("FINISH INIT PROCEDURE")
            return A, b, c_, N, B, v, msg
        else:
            msg = "No solution"
            return A, b, c, N, B, v, msg

    @staticmethod
    def _get_entering_var(c, N):
        var = None
        min_i = max(N) + 1
        for i, coeff in enumerate(c):
            if (coeff > 0) and (N[i] < min_i):
                var = N[i]
                min_i = N[i]
        return var

    @staticmethod
    def _get_leaving_var(A, b, N, B, e):
        e_i = N.index(e)

        var = None
        delta = []
        for i, row in enumerate(A):
            coeff = row[e_i]
            if coeff > 0:
                delta.append(b[i] / coeff)
            else:
                delta.append(float("inf"))

        min_i = max(B) + 1
        min_delta = float("inf")
        for i, d in enumerate(delta):
            if (d < float("inf")) and (d <= min_delta):
                if d < min_delta:
                    min_delta = d
                    min_i = B[i]
                else:
                    if B[i] < min_i:
                        min_delta = d
                        min_i = B[i]
        if min_delta < float("inf"):
            var = min_i

        return var

    @staticmethod
    def _pivot(A, b, c, v, N, B, l, e):
        """
        Pivot the system.

        l: leaving the basis variable index
        e: entering the basis variable index
        """
        leaving_row_i = B.index(l)
        entering_col_i = N.index(e)

        coeff = A[leaving_row_i][entering_col_i]
        A[leaving_row_i][entering_col_i] = 1
        A[leaving_row_i] = [x / coeff for x in A[leaving_row_i]]
        b[leaving_row_i] /= coeff

        for i in range(len(A)):
            if i == leaving_row_i:
                continue

            coeff = A[i][entering_col_i]
            b[i] = b[i] - coeff * b[leaving_row_i]
            for j in range(len(A[0])):
                A[i][j] = A[i][j] - coeff * A[leaving_row_i][j]
            A[i][entering_col_i] = -1 * coeff * A[leaving_row_i][entering_col_i]

        coeff = c[entering_col_i]
        v += coeff * b[leaving_row_i]
        for j in range(len(c)):
            c[j] = c[j] - coeff * A[leaving_row_i][j]
        c[entering_col_i] = -1 * coeff * A[leaving_row_i][entering_col_i]

        B[leaving_row_i] = e
        N[entering_col_i] = l

        return A, b, c, v, N, B


if __name__ == "__main__":
    n, m = list(map(int, input().split()))
    A = []
    for i in range(n):
        A += [list(map(int, input().split()))]
    b = list(map(int, input().split()))
    c = list(map(int, input().split()))

    aa = AdAllocation(A, b, c)
    aa.solve()
