#include <fstream>
#include <string>
#include <iostream>
#include <queue>  
#include <vector> 
#include <unordered_map>
#include <stack> 
#include <set> 
#include <algorithm>
#include <list>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct edge {
    int source;
    int target;
    int slength;
    int tlength;
    int overlap;
    char sorient;
    char torient;
}; // struct edge

void write_component_mappings(const std::unordered_map<int, int>& component_mapping, std::ofstream& outfile) {
  for(std::unordered_map<int, int>::const_iterator it = component_mapping.begin(); it != component_mapping.end(); ++it)
  {
      outfile << it->first << "\t" << it->second << std::endl;  
  }
}

inline bool containment(edge& curr_edge) {
  return (curr_edge.overlap >= curr_edge.tlength);
}

void write_edges(const std::vector<edge>& edges, std::ofstream& outfile) {
  for (edge curr_edge : edges) {
    outfile << curr_edge.source << "\t" << curr_edge.target << std::endl; 
  }
}

inline bool is_incoming_edge(const edge& curr_edge, int curr_node) {
  return (curr_edge.target == curr_node);
}

// This function that generates all pairs was largely taken from here: 
// https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
void combinationUtil(int arr[], int data[], int start, int end,
                     int index, int r, std::vector<std::pair<int,int> >& pairs)
{
    if (index == r)
    {
        std::pair<int,int> current_pair;
        current_pair = std::make_pair(data[0], data[1]);
        pairs.push_back(current_pair);
        return;
    }
 
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, r, pairs);
    }
}

std::vector<std::pair<int,int> > generate_all_pairs(std::vector<int> neighbors) {
  int combination_size = 2; 
  int pair[2]; 
  int num_neighbors = neighbors.size();
  std::vector<std::pair<int,int> > pairs;
  combinationUtil(neighbors.data(), pair, 0, num_neighbors-1, 0, combination_size, pairs); 
  
  return pairs;
}

void erase_edge(std::list<edge>& edges, int to_remove, bool source) {
  for (auto it = edges.begin(); it != edges.end(); it++) {
    edge curr_edge = *it;
    if (source) {
      if (curr_edge.source == to_remove) {
        edges.erase(it++);
      }
    }
    else {
      if (curr_edge.target == to_remove) {
        edges.erase(it++);
      }
    }
  }
}

bool anchored(const std::vector<int>&  pot_neighbors, 
  const std::vector<int>& curr_neighors, int curr_read, const std::set<int>& removed_nodes) {
  for (int neighbor : pot_neighbors) {
    if (neighbor == curr_read) continue; 
    if (removed_nodes.find(neighbor) != removed_nodes.end()) continue; 
    if (std::find(curr_neighors.begin(), curr_neighors.end(), neighbor) == curr_neighors.end()) {
      return false;
    }
  }
  return true; 
}

bool anchored_containment(const std::vector<int>&  pot_neighbors, 
  const std::vector<int>& curr_in_neighors, const std::vector<int>& curr_out_neighors, const std::vector<int>& possible_contained_nodes, int curr_read) {
  for (int neighbor : pot_neighbors) {
    if (neighbor == curr_read) continue; 
    bool not_in_in_neighbors = std::find(curr_in_neighors.begin(), curr_in_neighors.end(), neighbor) == curr_in_neighors.end();
    bool not_in_out_neighbors = std::find(curr_out_neighors.begin(), curr_out_neighors.end(), neighbor) == curr_out_neighors.end(); 
    bool not_in_possible_contained_nodes = std::find(possible_contained_nodes.begin(), possible_contained_nodes.end(), neighbor) == possible_contained_nodes.end(); 
    if (not_in_in_neighbors && not_in_out_neighbors && not_in_possible_contained_nodes) {
      return false;
    }
  }
  return true; 
}

bool possible_articulation_point_containment(int potential_contained_node, int curr_read, 
  const std::unordered_map<int, std::vector<int> >& trans_adjacency_list, 
  const std::unordered_map<int, std::vector<int> >& trans_comp, 
  const std::vector<int>& curr_incoming_nodes, const std::vector<int>& curr_outgoing_nodes, const std::vector<int>& possible_contained_nodes) 
{ 
  // Get all of the neighbors of the contained nodes
  std::vector<int> pot_contained_out_neighbors;
  std::vector<int> pot_contained_in_neighbors;
  bool has_outgoing_edges = trans_adjacency_list.find(potential_contained_node) != trans_adjacency_list.end();
  bool has_incoming_edges = trans_comp.find(potential_contained_node) != trans_comp.end();
  if (has_outgoing_edges) {
    pot_contained_out_neighbors = trans_adjacency_list.at(potential_contained_node); 
  }
  if (has_incoming_edges) {
    pot_contained_in_neighbors = trans_comp.at(potential_contained_node); 
  }

  bool incoming_neighbor_anchored = anchored_containment(pot_contained_in_neighbors, curr_incoming_nodes, curr_outgoing_nodes, possible_contained_nodes, curr_read); 
  bool outgoing_neighbor_anchored = anchored_containment(pot_contained_out_neighbors, curr_incoming_nodes, curr_outgoing_nodes, possible_contained_nodes, curr_read); 
  bool not_articulation_point = incoming_neighbor_anchored && outgoing_neighbor_anchored; 
  return not_articulation_point;
  
}


bool possible_articulation_point(int potential_trans_node, int curr_read, int partner_read,
  const std::unordered_map<int, std::vector<int> >& trans_adjacency_list, 
  const std::unordered_map<int, std::vector<int> >& trans_comp, 
  const std::vector<int>& curr_incoming_nodes, const std::vector<int>& curr_outgoing_nodes, const std::set<int>& removed_nodes) { 
  std::vector<int> pot_trans_out_neighbors;
  std::vector<int> pot_trans_in_neighbors;
  bool has_outgoing_edges = trans_adjacency_list.find(potential_trans_node) != trans_adjacency_list.end();
  bool has_incoming_edges = trans_comp.find(potential_trans_node) != trans_comp.end();
  if (has_outgoing_edges) {
    pot_trans_out_neighbors = trans_adjacency_list.at(potential_trans_node); 
  }
  if (has_incoming_edges) {
    pot_trans_in_neighbors = trans_comp.at(potential_trans_node); 
  }

  bool incoming_neighbor_anchored = (anchored(pot_trans_in_neighbors, curr_incoming_nodes, curr_read, removed_nodes) || anchored(pot_trans_in_neighbors, curr_incoming_nodes, partner_read, removed_nodes)); 
  bool outgoing_neighbor_anchored = (anchored(pot_trans_out_neighbors, curr_outgoing_nodes, curr_read, removed_nodes) || anchored(pot_trans_out_neighbors, curr_outgoing_nodes, partner_read, removed_nodes)); 
  
  return (incoming_neighbor_anchored && outgoing_neighbor_anchored); 
}

// This function will try to remove transitive nodes without disconnecting the graph
// This will return the number of edges it had to process
int find_non_articulate_transitive_nodes(std::vector<edge>& to_process,  int curr_read, 
  std::unordered_map<int, std::vector<int> >& trans_adjacency_list, std::unordered_map<int, std::vector<int> >& trans_comp, 
  std::unordered_map<int, int>&  component_mapping,
  std::set<int>& removed_nodes, std::list<edge>& curr_incoming_edges, std::list<edge>& curr_outgoing_edges) {

  // This classifies each edge as either an coming edge or an outgoing edge
  // Change this into a for loop. It is not necessary to delete the nodes in to_process. 
  // We just have to remember to clear it whenever we are returning. 
  std::vector<int> incoming_nodes; 
  std::vector<int> outgoing_nodes; 
  std::vector<int> possible_contained_nodes;
  while(!to_process.empty()) {
    edge curr_edge = to_process.back();
    if (containment(curr_edge)) {
      possible_contained_nodes.push_back(curr_edge.target);
      // If the current read is contained within another, we just remove it. It can't provide transitive closures 
      if (curr_edge.target == curr_read) {
        removed_nodes.clear();
        removed_nodes.insert(curr_read);
        to_process.clear();
        curr_incoming_edges.clear();
        curr_outgoing_edges.clear();
        return 0; 
      }
    }
    else {
      if (is_incoming_edge(curr_edge, curr_read)) {
        curr_incoming_edges.push_back(curr_edge);
        incoming_nodes.push_back(curr_edge.source);
      }
      else {
        curr_outgoing_edges.push_back(curr_edge);
        outgoing_nodes.push_back(curr_edge.target);
      }
    }
    to_process.pop_back();
  }
  
  for(int possible_contained_node : possible_contained_nodes) {
    if (possible_articulation_point_containment(possible_contained_node, curr_read, trans_adjacency_list, trans_comp, incoming_nodes, outgoing_nodes, possible_contained_nodes)) {
      removed_nodes.insert(possible_contained_node);
      erase_edge(curr_outgoing_edges, possible_contained_node, false); 
    }
    else {
      edge poss_contained_edge;
      poss_contained_edge.source = curr_read; 
      poss_contained_edge.target = possible_contained_node;
      curr_outgoing_edges.push_back(poss_contained_edge);
      outgoing_nodes.push_back(poss_contained_edge.target);
    }
  }

  int possible_edges_to_process = curr_incoming_edges.size() + curr_outgoing_edges.size();
  // See if current read forms an overlap that contains another read 
  int num_incoming = curr_incoming_edges.size();
  // if there are at least two incoming nodes, check to see if any of the node share an edge
  if (num_incoming >= 2) {
    std::vector<std::pair<int,int> > incoming_pairs = generate_all_pairs(incoming_nodes);
    for (std::pair<int, int> pair : incoming_pairs) {
      // Read1 and read2 have outgoing reads to the current read. Whichever has an incoming edge gets removed
      int read1 = pair.first;
      int read2 = pair.second; 
      // If either of them have been removed, then just proceed 
      if ((removed_nodes.find(read1) != removed_nodes.end()) || (removed_nodes.find(read2) != removed_nodes.end())) continue; 
      // Check if read1 has outgoing edges
      bool read1_outgoing = trans_adjacency_list.find(read1) != trans_adjacency_list.end(); 
      // Check if read2 has outgoing edges
      bool read2_outgoing = trans_adjacency_list.find(read2) != trans_adjacency_list.end(); 
      // Check if read1 has an edge to read2
      if (read1_outgoing && 
        (find(trans_adjacency_list.at(read1).begin(), trans_adjacency_list.at(read1).end(), read2) 
          != trans_adjacency_list.at(read1).end()) ) {
        if (possible_articulation_point(read2, curr_read, read1, trans_adjacency_list, trans_comp, incoming_nodes, outgoing_nodes, removed_nodes)) {
          removed_nodes.insert(read2); 
          erase_edge(curr_incoming_edges, read2, true); 
        }
      }
      else if (read2_outgoing && 
        (find(trans_adjacency_list.at(read2).begin(), trans_adjacency_list.at(read2).end(), read1) 
          != trans_adjacency_list.at(read2).end())) {
        if (possible_articulation_point(read1, curr_read, read2, trans_adjacency_list, trans_comp, incoming_nodes, outgoing_nodes, removed_nodes)) {
          removed_nodes.insert(read1); 
          erase_edge(curr_incoming_edges, read1, true); 
        }
      }
    }
  }
  
  int num_outgoing = curr_outgoing_edges.size();
  
  if (num_outgoing >= 2) {
    // Read1 and read2 have incoming edges from the current read. Whichever has an outgoing edge to another is the one that gets removed
    std::vector<std::pair<int,int> > outgoing_pairs = generate_all_pairs(outgoing_nodes);
    for (std::pair<int, int> pair : outgoing_pairs) {
      int read1 = pair.first;
      int read2 = pair.second; 
      // If either of them have been removed, then just proceed 
      if ((removed_nodes.find(read1) != removed_nodes.end()) || (removed_nodes.find(read2) != removed_nodes.end())) continue; 
      int read1_outgoing = trans_adjacency_list.find(read1) != trans_adjacency_list.end(); 
      int read2_outgoing = trans_adjacency_list.find(read2) != trans_adjacency_list.end(); 
      if (read1_outgoing && 
        (find(trans_adjacency_list.at(read1).begin(), trans_adjacency_list.at(read1).end(), read2) != trans_adjacency_list.at(read1).end())) {
        if (possible_articulation_point(read1, curr_read, read2, trans_adjacency_list, trans_comp, incoming_nodes, outgoing_nodes,removed_nodes)) {
          removed_nodes.insert(read1); 
          erase_edge(curr_outgoing_edges, read1, false); 
        }
 
      }
      else if (read2_outgoing && 
        (find(trans_adjacency_list.at(read2).begin(), trans_adjacency_list.at(read2).end(), read1) != trans_adjacency_list.at(read2).end())) { 
        if (possible_articulation_point(read2, curr_read, read1, trans_adjacency_list, trans_comp, incoming_nodes, outgoing_nodes,removed_nodes)) {
          removed_nodes.insert(read2); 
          erase_edge(curr_outgoing_edges, read2, false); 
        }
      }
    }
  }

  // If the new incoming read is joining multiple components, 
  // Don't consider removing it
  std::set<int> component_roots; 
  for(int node : incoming_nodes) {
    int root = component_mapping[node];
    if (component_roots.find(root) == component_roots.end()) {
      component_roots.insert(root);
      if(component_roots.size() > 1) {
        return possible_edges_to_process; 
      }
    }
  }

  for(int node : outgoing_nodes) {
    int root = component_mapping[node];
    if (component_roots.find(root) == component_roots.end()) {
      component_roots.insert(root);
      if(component_roots.size() > 1) {
        return possible_edges_to_process; 
      }
    }
  }

  if (removed_nodes.size() == 0) {
    num_incoming = curr_incoming_edges.size();
    num_outgoing = curr_outgoing_edges.size();
    // See if the current read is  contained with an overlap of two other reads
    if ((num_incoming > 0) && (num_outgoing > 0)) {
      for (edge incoming_edge : curr_incoming_edges) {
        for (edge outgoing_edge : curr_outgoing_edges) {
          // The source of the incoming edge is pointing at the current read
          // The current read is pointing to the target of the outgoing edge
          // If the source of the incoming read shares an edge with the target of the outgoing edge, remove the current node
          int read1 = incoming_edge.source; 
          int read2 = outgoing_edge.target; 
          bool read1_has_edge_to_read2 = find(trans_adjacency_list.at(read1).begin(), trans_adjacency_list.at(read1).end(), read2) != trans_adjacency_list.at(read1).end();
          if (read1_has_edge_to_read2) {
            removed_nodes.insert(curr_read);
            curr_incoming_edges.clear();
            curr_outgoing_edges.clear();
            return possible_edges_to_process; 
          }
        }
      }
    }
  }
  return possible_edges_to_process;
}

void remove_from_graph(std::unordered_map<int, std::vector<int> >& graph, 
                  std::unordered_map<int, std::vector<int> >& graph_complement, int to_remove) {

  // Using one graph to update the other graph
  // For example, if graph were the regular adjacency list and graph complement were the inverse, 
  // Then when we remove the read r1, we would have to erase r1 from the original graph and use the 
  // Graph complement to update any reads that have an outgoing edge to r1 
  std::vector<int> curr_edges = graph[to_remove]; 
  for (int read : curr_edges) {
    if (graph_complement.find(read) != graph_complement.end()) {
      //std::vector<int> curr_comp_edges = graph_complement[read];
      graph_complement[read].erase(std::remove(graph_complement[read].begin(), graph_complement[read].end(), to_remove), graph_complement[read].end());
    }
  }
  graph.erase(to_remove);

}

void update_edges(const std::set<int>& removed_nodes, 
                  std::unordered_map<int, std::vector<int> >& trans_comp, 
                  std::unordered_map<int, std::vector<int> >& trans_adjacency_list) {

  for (std::set<int>::iterator it = removed_nodes.begin(); it != removed_nodes.end(); ++it) {
    int current_removed_node = *it;
    // Remove all of our outgoing edges
    if (trans_adjacency_list.find(current_removed_node) != trans_adjacency_list.end()) {
      remove_from_graph(trans_adjacency_list, trans_comp, current_removed_node);
    }
    // Remove all incoming edges 
    if (trans_comp.find(current_removed_node) != trans_comp.end()) {
      remove_from_graph(trans_comp, trans_adjacency_list, current_removed_node);
    }
  }
}
inline void update_complete_graph(std::unordered_map<int, std::vector<int> >& complete_adjacency_list, 
                            std::set<int>& all_nodes, edge curr_edge) {
  
  // TODO: Do you want only nodes that have edges to be considered in all_nodes? 
  // Currently nodes that do not have edges are not included.
  all_nodes.insert(curr_edge.source);
  all_nodes.insert(curr_edge.target);
  
  // The source node has never been seen before 
  if (complete_adjacency_list.find(curr_edge.source) == complete_adjacency_list.end()) {
    complete_adjacency_list[curr_edge.source] = std::vector<int>();
  }
  complete_adjacency_list[curr_edge.source].push_back(curr_edge.target);

  // Comment the below lines out if you want complete_adjacency_list to be directed
  // The target node has never been seen before 
  if (complete_adjacency_list.find(curr_edge.target) == complete_adjacency_list.end()) {
    complete_adjacency_list[curr_edge.target] = std::vector<int>();
  }
  complete_adjacency_list[curr_edge.target].push_back(curr_edge.source);
}

inline void update_nodes(std::set<int>& irreducible_nodes, const std::set<int>& new_removed_nodes) {
  for (int removed_node : new_removed_nodes) {
    irreducible_nodes.erase(removed_node);
  }
}

void relabel_components(int curr_read, const bool incoming, const std::list<edge>& edge_list, 
                        std::unordered_map<int, int>&  mapping_components, 
                        std::unordered_map<int, std::set<int> >&  component_members, 
                        const std::set<int>& removed_nodes) {

  // Find the minimum root id
  int min_root = mapping_components[curr_read]; 
  for(edge curr_edge : edge_list) {
    int neighbor = curr_edge.source; 
    if (!incoming) neighbor = curr_edge.target;
    int neighbor_root = mapping_components[neighbor];
    if (removed_nodes.find(neighbor_root) != removed_nodes.end()) continue;  
    if (neighbor_root < min_root) min_root = neighbor_root;
  }

  // Handle the current read 
  component_members[min_root].insert(curr_read);
  if(min_root != mapping_components[curr_read]) {
    int prev_root = mapping_components[curr_read];
    for(int member : component_members[prev_root]) {
      mapping_components[member] = min_root;
      component_members[min_root].insert(member);
    }
    component_members.erase(prev_root);
  }
  
  // Handle all of the current read's neighbors 
  for(edge curr_edge : edge_list) {
    int neighbor = curr_edge.source; 
    if (!incoming) neighbor = curr_edge.target; 
    if (mapping_components[neighbor] == min_root) continue; 
    int prev_root = mapping_components[neighbor];
    for(int member : component_members[prev_root]) {
      mapping_components[member] = min_root;
      component_members[min_root].insert(member);
    }
    component_members.erase(prev_root);
  }
}

void update_components(int curr_read, const std::list<edge>& incoming_edges, const std::list<edge>& outgoing_edges, 
                        const std::set<int>& removed_nodes, std::unordered_map<int, int>&  mapping_components, 
                        std::unordered_map<int, std::set<int> >&  component_members, 
                         std::unordered_map<int, std::vector<int> >& trans_adjacency_list) {

  if (removed_nodes.find(curr_read) != removed_nodes.end()) {
    return;
  } 
  
  // The current read starts off as its own component 
  mapping_components[curr_read] = curr_read; 
  component_members[curr_read] = std::set<int>(); 
  component_members[curr_read].insert(curr_read);

  relabel_components(curr_read, true, incoming_edges, mapping_components, component_members, removed_nodes);
  relabel_components(curr_read, false, outgoing_edges, mapping_components, component_members, removed_nodes);


  // Updates nodes affected by transitive closures
  for(int removed_node : removed_nodes) {
    int root_of_removed_node = mapping_components[removed_node]; 
    if (component_members.find(root_of_removed_node) != component_members.end()) {
      component_members[root_of_removed_node].erase(removed_node);
    }

    if (component_members.find(removed_node) != component_members.end()) {
      for (auto removed_node_member : component_members[removed_node]) {
        if (removed_nodes.find(removed_node_member) == removed_nodes.end()) {
          component_members[mapping_components[curr_read]].insert(removed_node_member);
          mapping_components[removed_node_member] = mapping_components[curr_read];
        }
      }
      component_members.erase(removed_node);
    }

    mapping_components.erase(removed_node); 
  }
}

void degree_stats(std::unordered_map<int, std::vector<int> >& irreducible_edges, std::unordered_map<int, std::vector<int>>& trans_comp, std::set<int>& irreducible_nodes, const std::string& output_file_prefix) {
  std::ofstream degree_file;
  std::string degree_file_name = output_file_prefix + ".irreducible_degree";
  degree_file.open (degree_file_name); 
  for (auto it = irreducible_nodes.begin(); it != irreducible_nodes.end(); ++it) {
    int node = *it; 
    int in_degree = 0;
    int out_degree = 0;
    if (trans_comp.find(node) != trans_comp.end()) {
      in_degree = trans_comp[node].size();
    }
    if (irreducible_edges.find(node) != irreducible_edges.end()) {
      out_degree = irreducible_edges[node].size();
    }
    degree_file << node << '\t' << in_degree << '\t' << out_degree << std::endl; 
  }
  degree_file.close();
}

int trans_closures_stream(std::unordered_map<int, std::vector<int> >& trans_adjacency_list, 
                      std::unordered_map<int, std::vector<int> >& complete_adjacency_list, 
                        std::set<int>& irreducible_nodes, std::set<int>& all_nodes, std::set<int>& removed_nodes,
                        std::unordered_map<int, int>&  component_mapping, std::unordered_map<int, std::set<int> >& component_members,
                        const std::string& output_file_prefix) {

  std::string storage_rate_filename = output_file_prefix + ".storage_rate";
  std::string component_convergence_filename = output_file_prefix + ".component_convergence";
  std::string removal_rate_filename = output_file_prefix + ".removal_rate";

  std::ofstream component_file;
  component_file.open (component_convergence_filename); 

  std::ofstream storage_file;
  storage_file.open (storage_rate_filename);

  std::ofstream removal_rate_file;
  removal_rate_file.open (removal_rate_filename);

  std::vector<edge> to_process;
  edge curr_edge; 
  // Transitive adjacency list complement 
  std::unordered_map<int, std::vector<int> > trans_comp; 
  int curr_read = -1;
  //This while condition loop is setting the variable values 
  int num_reads = 0; 
  int max_number_of_reads = 0; 

  while(std::cin >> curr_edge.source) {
    if (curr_edge.source == -1) {
      std::cin >> curr_read;
      all_nodes.insert(curr_read);
      irreducible_nodes.insert(curr_read);
      std::set<int> new_removed_nodes; 
      std::list<edge> incoming_edges; 
      std::list<edge> outgoing_edges; 

      int num_edges_to_process = to_process.size();
      // This line documents the number of components there are in the graph
      component_file << component_members.size() << std::endl; 

      int possible_edges_processed = find_non_articulate_transitive_nodes(to_process, curr_read, trans_adjacency_list, trans_comp, component_mapping, new_removed_nodes, incoming_edges, outgoing_edges);

      int num_edges_after_processing = incoming_edges.size() + outgoing_edges.size();
      // Collect stats on how many edges were removed for a node
      // removal_rate_file << curr_read << '\t' << complete_adjacency_list[curr_read].size() << '\t' << num_edges_to_process << '\t' << num_edges_after_processing << std::endl; 
      
      // There is a distinction between this line and the line prior because contained nodes count towards the degree
      removal_rate_file << curr_read << '\t' << complete_adjacency_list[curr_read].size() << '\t' << possible_edges_processed << '\t' << num_edges_after_processing << std::endl; 
      
      // After find containment, all that is left in incoming_edges and outgoing edges are the nodes that were added to the graph. Not the ones to be deleted 
      // Update the component mapping before we remove the edges of removed nodes and update the graph
      update_components(curr_read, incoming_edges, outgoing_edges, new_removed_nodes, component_mapping, component_members,trans_adjacency_list); 
      // Incrementally update the graph. We can change this so that it updates it in batches 
      update_edges(new_removed_nodes, trans_comp, trans_adjacency_list); 
      update_nodes(irreducible_nodes, new_removed_nodes);
      removed_nodes.insert(new_removed_nodes.begin(), new_removed_nodes.end()); 
      curr_read = -1; 
      ++num_reads;
      
      int num_reads_stored = irreducible_nodes.size();

      // This line documents the number of reads stored 
      storage_file << num_reads_stored << std::endl; 
      
      //The following line was used to get the maximum number of reads stored at any given point (ie irreducible_nodes)
      if (num_reads_stored > max_number_of_reads) max_number_of_reads = irreducible_nodes.size();
    }
    else {
      std::cin >> curr_edge.target >> curr_edge.slength >> curr_edge.tlength >> curr_edge.overlap >> curr_edge.sorient >> curr_edge.torient;
      // Forming complete graph 
      // Uncomment this if you want information on the complete graph 
      if (curr_edge.source == curr_edge.target) continue; 
      update_complete_graph(complete_adjacency_list, all_nodes, curr_edge);
      bool source_removed = removed_nodes.find(curr_edge.source) != removed_nodes.end();
      bool target_removed = removed_nodes.find(curr_edge.target) != removed_nodes.end();
      if (source_removed || target_removed) {
        // Ignore the edges that contain nodes that we have removed for transitive closures 
        continue; 
      } 
      irreducible_nodes.insert(curr_edge.source);
      irreducible_nodes.insert(curr_edge.target);
      bool source_has_no_outgoing_edges = trans_adjacency_list.find(curr_edge.source) == trans_adjacency_list.end(); 
      bool target_has_no_incoming_edges = trans_comp.find(curr_edge.target) == trans_comp.end();
      if (source_has_no_outgoing_edges) {
        trans_adjacency_list[curr_edge.source] = std::vector<int>();
      }

      if (target_has_no_incoming_edges) {
        trans_comp[curr_edge.target] = std::vector<int>();
      }   
  
      trans_adjacency_list[curr_edge.source].push_back(curr_edge.target);
      trans_comp[curr_edge.target].push_back(curr_edge.source); 

      to_process.push_back(curr_edge); 
    }
    
  }

  component_file.close();
  storage_file.close();
  degree_stats(trans_adjacency_list, trans_comp, irreducible_nodes, output_file_prefix);
  return max_number_of_reads; 
}

std::vector<edge> list_edges(std::unordered_map<int, std::vector<int> >& adjacency_list) {
  std::vector<edge> edges;
  for(std::unordered_map<int, std::vector<int> >::const_iterator it = adjacency_list.begin(); it != adjacency_list.end(); ++it)
  {
      int parent = it->first; 
      std::vector<int> children = it->second; 
      for (int child : children) {
        edge curr_edge;
        curr_edge.source = parent; 
        curr_edge.target = child;
        edges.push_back(curr_edge);
      }
  }
  return edges;

}

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options");
  // First parameter describes option name/short name
  // The second is parameter to option
  // The third is description
  desc.add_options()
  ("help,h", "print usage message")
  ("output,o", po::value<std::string>(), "Output file names")
  ("component_size,c", po::value<int>(), "Component Threshold")
  ("irreducible,i", po::value<std::string>(), "Irreducible Edge file name to write to")
  ("threshold,t", po::value<int>(), "Overlap Threshold");
  
  po::variables_map vm;        
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {  
      std::cout << desc << "\n";
      return 0;
  }
  std::string output_file_prefix = "default_out";
  if (vm.count("output")) {
      output_file_prefix = vm["output"].as<std::string>(); 
      std::cout << "Output file is " 
           << output_file_prefix << "." << std::endl;
  } else {
      std::cout << "Output file was not passed in, default is: " << output_file_prefix << std::endl;
  }

  int component_size_threshold = 10; 
  if (vm.count("component_size")) {
      component_size_threshold = vm["component_size"].as<int>(); 
      std::cout << "Component size threshold is set as " 
           << component_size_threshold << "." << std::endl;
  } else {
      std::cout << "Component size threshold was not passed in, default is: " << component_size_threshold << std::endl;
  }

  std::string stat_filename = output_file_prefix + ".trans_stats";
  std::ofstream stat_file;
  stat_file.open (stat_filename);

  std::string componet_size_filename = output_file_prefix + ".component_sizes";
  std::ofstream componet_size_file;
  componet_size_file.open (componet_size_filename);
  
  // This adjacency list will have transitive closures performed on it 
  std::unordered_map<int, std::vector<int> > trans_adjacency_list; 
  
  // This adjacency will not have any transitive closures done. This is an undirected graph. 
  // Change update_complete_graph if you want to make it directed 
  std::unordered_map<int, std::vector<int> > complete_adjacency_list; 
  std::set<int> irreducible_nodes; 
  std::set<int> all_nodes; 
  std::set<int> removed_nodes; 
  std::unordered_map<int, int> component_mapping;
  std::unordered_map<int, std::set<int> > component_members;  
  
  // The below call does a majority of the work 
  int max_number_of_reads = trans_closures_stream(trans_adjacency_list, complete_adjacency_list, irreducible_nodes, all_nodes,removed_nodes, component_mapping,component_members, output_file_prefix);

  /*
    The below lines just print out stats. Processing is already finished.
  */

  stat_file << "reads processed: " << all_nodes.size() << std::endl;
  stat_file << "Number of nodes removed: " << removed_nodes.size() << std::endl;
  stat_file << "Number of mapped nodes " << component_mapping.size() <<std::endl; 
  stat_file << "Number of trans nodes " << irreducible_nodes.size() << std::endl; 
  stat_file << "Number of connected components" << std::endl; 
  stat_file << component_members.size() << std::endl; 
  int num_mapped_nodes = 0; 
  int number_of_significant_components = 0; 
  std::vector< std::pair<int, int> > significant_components;
  std::set<int> nodes_with_membership; 
  stat_file << "Significant components " << std::endl; 
  for (auto member_pair : component_members) {
    int member = member_pair.first; 
    int component_size = member_pair.second.size();
    num_mapped_nodes+=component_size;
    nodes_with_membership.insert(member_pair.second.begin(), member_pair.second.end());
    if (component_size > component_size_threshold) {
      stat_file << member << " : " << component_size << std::endl;  
      ++number_of_significant_components;
    }
    componet_size_file << member << "\t" << component_size << std::endl;
  }
  componet_size_file.close();
  stat_file << "Number of significant components" << std::endl; 
  stat_file << number_of_significant_components << std::endl;
  stat_file << "Number of component members" << std::endl; 
  stat_file << num_mapped_nodes << std::endl; 
  stat_file << "Maximum number of nodes stored" << std::endl; 
  stat_file << max_number_of_reads << std::endl; 

  stat_file.close();

  std::string irreducible_edge_filename = output_file_prefix + ".irreducible_edges";
  std::ofstream irreducible_file;
  irreducible_file.open(irreducible_edge_filename); 
  std::vector<edge> edges = list_edges(trans_adjacency_list);
  write_edges(edges, irreducible_file);
  irreducible_file.close();

  std::string trans_mapping_filename = output_file_prefix + ".trans_mapping";
  std::ofstream trans_mapping_file;
  trans_mapping_file.open(trans_mapping_filename); 
  write_component_mappings(component_mapping, trans_mapping_file);
  trans_mapping_file.close();

  std::cout << "Using Boost "     
          << BOOST_VERSION / 100000     << "."  // major version
          << BOOST_VERSION / 100 % 1000 << "."  // minor version
          << BOOST_VERSION % 100                // patch level
          << std::endl;
  
  return 0;
} // main