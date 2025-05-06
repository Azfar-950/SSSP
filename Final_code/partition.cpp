#include <metis.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <stdexcept>

void read_graph(const std::string& filename, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, std::vector<idx_t>& adjwgt, std::vector<std::tuple<idx_t, idx_t, idx_t>>& edges) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open graph file: " + filename);
    }

    std::vector<std::vector<std::pair<idx_t, idx_t>>> adj_list;
    idx_t u, v, w;
    idx_t max_vertex = -1;

    while (file >> u >> v >> w) {
        if (u < 0 || v < 0) {
            throw std::runtime_error("Negative vertex ID found");
        }
        if (w < 0) {
            throw std::runtime_error("Negative weight found at edge (" + std::to_string(u) + ", " + std::to_string(v) + ")");
        }
        max_vertex = std::max(max_vertex, std::max(u, v));
        if (adj_list.size() <= static_cast<size_t>(max_vertex)) {
            adj_list.resize(max_vertex + 1);
        }
        adj_list[u].emplace_back(v, w);
        adj_list[v].emplace_back(u, w);
        edges.emplace_back(u, v, w);
    }

    idx_t nvtxs = adj_list.size();
    xadj.push_back(0);
    for (idx_t i = 0; i < nvtxs; ++i) {
        for (const auto& [v, w] : adj_list[i]) {
            adjncy.push_back(v);
            adjwgt.push_back(w);
        }
        xadj.push_back(adjncy.size());
    }

    std::cout << "Read graph with " << nvtxs << " vertices and " << adjncy.size() / 2 << " edges\n";
}

void extract_subgraphs(const std::vector<idx_t>& part, const std::vector<std::tuple<idx_t, idx_t, idx_t>>& edges, idx_t nparts) {
    std::vector<std::vector<std::tuple<idx_t, idx_t, idx_t>>> subgraphs(nparts);
    std::vector<std::set<idx_t>> local_nodes(nparts);
    std::vector<std::set<idx_t>> ghost_nodes(nparts);

    for (const auto& [u, v, w] : edges) {
        idx_t u_part = part[u];
        idx_t v_part = part[v];
        if (u_part == v_part) {
            subgraphs[u_part].emplace_back(u, v, w);
            local_nodes[u_part].insert(u);
            local_nodes[u_part].insert(v);
        } else {
            subgraphs[u_part].emplace_back(u, v, w);
            local_nodes[u_part].insert(u);
            ghost_nodes[u_part].insert(v);
            subgraphs[v_part].emplace_back(u, v, w);
            local_nodes[v_part].insert(v);
            ghost_nodes[v_part].insert(u);
        }
    }

    for (idx_t part = 0; part < nparts; ++part) {
        std::ofstream out("subgraph_" + std::to_string(part) + ".txt");
        for (const auto& [u, v, w] : subgraphs[part]) {
            out << u << " " << v << " " << w << "\n";
        }
        out.close();
        std::cout << "Wrote subgraph_" << part << ".txt with " << subgraphs[part].size() << " edges\n";

        std::ofstream nodes_out("subgraph_" + std::to_string(part) + "_nodes.txt");
        nodes_out << "# Local nodes\n";
        for (idx_t node : local_nodes[part]) {
            nodes_out << node << "\n";
        }
        nodes_out << "# Ghost nodes\n";
        for (idx_t node : ghost_nodes[part]) {
            nodes_out << node << "\n";
        }
        nodes_out.close();
        std::cout << "Wrote subgraph_" << part << "_nodes.txt with "
                  << local_nodes[part].size() << " local nodes and "
                  << ghost_nodes[part].size() << " ghost nodes\n";
    }
}

int main() {
    // Graph data
    std::vector<idx_t> xadj, adjncy, adjwgt;
    std::vector<std::tuple<idx_t, idx_t, idx_t>> edges;
    try {
        read_graph("graph.txt", xadj, adjncy, adjwgt, edges);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    idx_t nvtxs = xadj.size() - 1;
    idx_t ncon = 1;
    idx_t nparts = 4; // Adjust based on desired partitions
    std::vector<idx_t> part(nvtxs);
    idx_t objval;

    // Validate edge weights
    for (idx_t w : adjwgt) {
        if (w < 0) {
            std::cerr << "Error: Negative edge weight detected in adjwgt: " << w << std::endl;
            return 1;
        }
    }

    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_UFACTOR] = 100;
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM;
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
    options[METIS_OPTION_NCUTS] = 1;
    options[METIS_OPTION_NITER] = 10;
    options[METIS_OPTION_DBGLVL] = 0;

    // Partition graph
    int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), nullptr, nullptr, adjwgt.data(),
                                  &nparts, nullptr, nullptr, options, &objval, part.data());
    if (ret != METIS_OK) {
        std::cerr << "METIS partitioning failed with error code: " << ret << "\n";
        return 1;
    }
    if (objval < 0) {
        std::cerr << "Error: METIS returned negative edge cut: " << objval << "\n";
        return 1;
    }

    // Output partitions
    std::cout << "Edge cut: " << objval << "\n";
    for (idx_t i = 0; i < nvtxs; ++i) {
        std::cout << "Vertex " << i << " -> Partition " << part[i] << "\n";
    }

    // Extract subgraphs
    extract_subgraphs(part, edges, nparts);

    return 0;
}