extern "C" __global__ void initialize_parents(unsigned int* parents, unsigned int n) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        parents[idx] = idx;
    }
}

extern "C" __global__ void union_step(unsigned int* parents, const uint2* edges, int num_edges) {
    int edge_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (edge_idx >= num_edges) return;

    uint2 edge = edges[edge_idx];
    unsigned int u = edge.x;
    unsigned int v = edge.y;

    while (true) {
        unsigned int root_u = u;
        while (parents[root_u] != root_u) {
            unsigned int parent = parents[root_u];
            parents[root_u] = parents[parent];
            root_u = parent;
        }

        unsigned int root_v = v;
        while (parents[root_v] != root_v) {
            unsigned int parent = parents[root_v];
            parents[root_v] = parents[parent];
            root_v = parent;
        }

        if (root_u == root_v) break;

        // Link smaller to larger to keep trees proper.
        if (root_u > root_v) {
            unsigned int tmp = root_u; root_u = root_v; root_v = tmp;
        }
        unsigned int old = atomicCAS(&parents[root_u], root_u, root_v);
        if (old == root_u) break;
        u = old;
        v = root_v;
    }
}

extern "C" __global__ void pointer_jump(unsigned int* parents, unsigned int n, int* changed) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        unsigned int p  = parents[idx];
        unsigned int gp = parents[p];
        if (p != gp) {
            parents[idx] = gp;
            atomicOr(changed, 1); // plain store would be a data race
        }
    }
}
