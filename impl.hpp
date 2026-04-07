// Implementation file injected by build

// matrix methods
matrix::matrix(int m_, int n_) { 
    if (m_ <= 0 || n_ <= 0) { m = n = 0; data = nullptr; return; }
    m = m_; n = n_;
    data = new fraction*[m];
    for (int i = 0; i < m; ++i) {
        data[i] = new fraction[n];
        for (int j = 0; j < n; ++j) data[i][j] = fraction(0);
    }
}

matrix::matrix(const matrix &obj) {
    if (obj.m <= 0 || obj.n <= 0 || obj.data == nullptr) { m = n = 0; data = nullptr; return; }
    m = obj.m; n = obj.n;
    data = new fraction*[m];
    for (int i = 0; i < m; ++i) {
        data[i] = new fraction[n];
        for (int j = 0; j < n; ++j) data[i][j] = obj.data[i][j];
    }
}

matrix::matrix(matrix &&obj) noexcept { m = obj.m; n = obj.n; data = obj.data; obj.m = obj.n = 0; obj.data = nullptr; }

matrix::~matrix() {
    if (data) { for (int i = 0; i < m; ++i) delete [] data[i]; delete [] data; }
    data = nullptr; m = n = 0;
}

matrix &matrix::operator=(const matrix &obj) {
    if (this == &obj) return *this;
    if (data) { for (int i = 0; i < m; ++i) delete [] data[i]; delete [] data; }
    data = nullptr; m = n = 0;
    if (obj.m <= 0 || obj.n <= 0 || obj.data == nullptr) return *this;
    m = obj.m; n = obj.n;
    data = new fraction*[m];
    for (int i = 0; i < m; ++i) {
        data[i] = new fraction[n];
        for (int j = 0; j < n; ++j) data[i][j] = obj.data[i][j];
    }
    return *this;
}

fraction &matrix::operator()(int i, int j) {
    if (i <= 0 || i > m || j < 0 || j >= n || data == nullptr) throw matrix_error();
    return data[i - 1][j];
}

matrix operator*(const matrix &lhs, const matrix &rhs) {
    if (lhs.n == 0 || rhs.m == 0 || lhs.m == 0 || rhs.n == 0 || lhs.n != rhs.m) throw matrix_error();
    matrix res(lhs.m, rhs.n);
    for (int i = 0; i < lhs.m; ++i) {
        for (int k = 0; k < lhs.n; ++k) {
            const fraction &lik = lhs.data[i][k];
            if (lik == fraction(0)) continue;
            for (int j = 0; j < rhs.n; ++j) {
                res.data[i][j] = res.data[i][j] + lik * rhs.data[k][j];
            }
        }
    }
    return res;
}

matrix matrix::transposition() {
    if (m == 0 || n == 0 || data == nullptr) throw matrix_error();
    matrix t(n, m);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j) t.data[j][i] = data[i][j];
    return t;
}

fraction matrix::determination() {
    if (m == 0 || n == 0 || data == nullptr || m != n) throw matrix_error();
    int size = n;
    std::vector<std::vector<fraction>> a(size, std::vector<fraction>(size));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) a[i][j] = data[i][j];
    bool sign_pos = true;
    for (int col = 0, row = 0; col < size && row < size; ++col, ++row) {
        int pivot = row;
        while (pivot < size && a[pivot][col] == fraction(0)) ++pivot;
        if (pivot == size) return fraction(0);
        if (pivot != row) { std::swap(a[pivot], a[row]); sign_pos = !sign_pos; }
        for (int i = row + 1; i < size; ++i) {
            if (a[i][col] == fraction(0)) continue;
            fraction factor = a[i][col] / a[row][col];
            for (int j = col; j < size; ++j) a[i][j] = a[i][j] - factor * a[row][j];
        }
    }
    fraction det(1);
    for (int i = 0; i < size; ++i) det = det * a[i][i];
    if (!sign_pos) det = fraction(0) - det;
    return det;
}

// resistive_network methods
resistive_network::resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
    : interface_size(interface_size_), connection_size(connection_size_),
      adjacency(connection_size_, interface_size_), conduction(connection_size_, connection_size_) {
    for (int e = 0; e < connection_size; ++e) {
        int u = from[e] - 1;
        int v = to[e] - 1;
        adjacency(e + 1, u) = adjacency(e + 1, u) + fraction(1);
        adjacency(e + 1, v) = adjacency(e + 1, v) + fraction(-1);
        fraction g = fraction(1) / resistance[e];
        conduction(e + 1, e) = g;
    }
}

static matrix build_laplacian(matrix &adjacency, matrix &conduction) {
    matrix At = adjacency.transposition();
    matrix CA = conduction * adjacency;
    return At * CA;
}

static matrix reduced_top_left(const matrix &L) {
    int n = const_cast<matrix&>(L).rows();
    matrix R(n - 1, n - 1);
    for (int i = 0; i < n - 1; ++i)
        for (int j = 0; j < n - 1; ++j)
            R(i + 1, j) = const_cast<matrix&>(L)(i + 1, j);
    return R;
}

static fraction cramer_component(matrix A, const std::vector<fraction> &b, int j_index) {
    fraction detA = A.determination();
    if (detA == fraction(0)) throw resistive_network_error();
    matrix Aj = A;
    int n = A.rows();
    for (int i = 0; i < n; ++i) Aj(i + 1, j_index) = b[i];
    fraction detAj = Aj.determination();
    return detAj / detA;
}

fraction resistive_network::get_equivalent_resistance(int interface_id1, int interface_id2) {
    matrix L = build_laplacian(adjacency, conduction);
    matrix Lr = reduced_top_left(L);
    int n = interface_size;
    std::vector<fraction> I1(n - 1, fraction(0));
    auto addI = [&](int idx, const fraction &val){ if (idx >= 0 && idx < n - 1) I1[idx] = I1[idx] + val; };
    addI(interface_id1 - 1, fraction(1));
    addI(interface_id2 - 1, fraction(-1));
    auto voltage_of = [&](int id)->fraction{
        if (id == n) return fraction(0);
        return cramer_component(Lr, I1, id - 1);
    };
    fraction vi = voltage_of(interface_id1);
    fraction vj = voltage_of(interface_id2);
    return vi - vj;
}

fraction resistive_network::get_voltage(int id, fraction current[]) {
    matrix L = build_laplacian(adjacency, conduction);
    matrix Lr = reduced_top_left(L);
    int n = interface_size;
    std::vector<fraction> I1(n - 1);
    for (int i = 0; i < n - 1; ++i) I1[i] = current[i];
    return cramer_component(Lr, I1, id - 1);
}

fraction resistive_network::get_power(fraction voltage[]) {
    fraction total(0);
    for (int e = 0; e < connection_size; ++e) {
        int u = -1, v = -1;
        for (int j = 0; j < interface_size; ++j) {
            fraction val = adjacency(e + 1, j);
            if (val == fraction(1)) u = j; else if (val == fraction(-1)) v = j;
        }
        fraction dV = voltage[u] - voltage[v];
        fraction g = conduction(e + 1, e);
        total = total + g * dV * dV;
    }
    return total;
}
