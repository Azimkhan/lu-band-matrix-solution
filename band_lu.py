from fractions import Fraction

def read_matrix(filename):
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        b = int(f.readline().strip())
        matrix = []
        for i in range(n):
            line = f.readline().strip()
            raw_nums = line.split()
            if len(raw_nums) != n:
                raise ValueError('row %d length must be %d numbers' % (i+1, n))
            row = tuple(Fraction.from_float(float(x)) for x in raw_nums)
            matrix.append(row)
        line = f.readline().strip()
        raw_nums = line.split()
        if len(raw_nums) != n:
            raise ValueError('c vector length must be %d numbers' % n)
        c = tuple(Fraction.from_float(float(x)) for x in raw_nums)

        return n,b,matrix,c


def lu_decomposition(m,b,n):
    # make list of lists from list of tuples
    u = list(list(x) for x in m)
    l = list(list(Fraction(1) if y == x else Fraction() for y in range(n)) for x in range(n))
    for k in range(n-1):
        for i in range(k+1, min([k+b, n])+1):
            l[i][k] = u[i][k] / u[k][k]
            for j in range(k, min([k+b,n])+1):
                u[i][j] = u[i][j] - l[i][k] * u[k][j]
    return l,u

def solve_lu(l,u,c,n):
    y = []
    for i in range(n):
        if i == 0:
            y.append(c[i])
        else:
            s = Fraction()
            for j in range(i):
                s += l[i][j]*y[j]
            y.append(c[i] - s)

    x = []
    for i in reversed(range(n)):
        if i == n-1:
            x.append(y[i]/u[i][i])
        else:
            s = Fraction()
            for j in reversed(range(i+1,n)):
                s += u[i][j] * x[n-j-1]
            x.append((y[i]-s)/u[i][i])
    return tuple(reversed(x))

# prints matrix
def print_2d(m, label=None, new_line=False):
    label_printed = False
    for row in m:
        row_content = print_1d(row,ret=True).strip()
        if label:
            if not label_printed:
                print('%s = %s' % (label, row_content))
                label_printed = True
            else:
                print('%s   %s' % (' ' * len(label), row_content)) 
        else:
            print(row_content)

    if new_line:
        print()

# prints vector
def print_1d(row,label=None, new_line=False, ret=False):
    s = []
    if label:
        s.append('%s = ' % label)
    s.append("| %s |\n" % "  ".join('%2.3f' % x for x in row))
    if new_line:
        s.append("\n")
    content = "".join(s)
    if ret:
        return content
    else:
        print(content,end='')

'''

This application sovles Ax = C linear equation system where
  A is a banded matrix
  C is the right side of the system
  n - matrix size
  b - half band size

The matrix file should look like:
  n
  b
  a11 a12 a13 ... a1n
  a21 a22 a22 ... a2n
  .
  .
  .
  an1 an2 an3 ... ann

'''
if __name__ == '__main__':
    import sys
    filename = sys.argv[1] if len(sys.argv) > 1 else 'matrix.txt'
    n,b,a,c = read_matrix(filename)
    print_2d(a,'A', True)
    print_1d(c,'C', True)

    l,u = lu_decomposition(a,b,n)

    x = solve_lu(l,u,c,n)
    print_1d(x, 'X',True)