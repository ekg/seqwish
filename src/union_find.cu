extern "C" __global__ void initialize_parents(int* parents, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        parents[idx] = idx;
    }
}

extern "C" __global__ void union_step(int* parents, const int2* edges, int num_edges) {
    int edge_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (edge_idx < num_edges) {
        int2 edge = edges[edge_idx];
        int u = edge.x;
        int v = edge.y;
        
        while (true) {
            int root_u = u;
            while (parents[root_u] != root_u) {
                int parent = parents[root_u];
                parents[root_u] = parents[parent];
                root_u = parent;
            }
            
            int root_v = v;
            while (parents[root_v] != root_v) {
                int parent = parents[root_v];
                parents[root_v] = parents[parent];
                root_v = parent;
            }
            
            if (root_u == root_v) {
                break;
            }
            
            if (root_u < root_v) {
                int old = atomicCAS(&parents[root_u], root_u, root_v);
                if (old == root_u) {
                    break;
                }
                u = old;
                v = root_v;
            } else {
                int old = atomicCAS(&parents[root_v], root_v, root_u);
                if (old == root_v) {
                    break;
                }
                u = root_u;
                v = old;
            }
        }
    }
}

extern "C" __global__ void pointer_jump(int* parents, int n, int* changed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        int p = parents[idx];
        int gp = parents[p];
        if (p != gp) {
            parents[idx] = gp;
            *changed = 1;
        }
    }
}
