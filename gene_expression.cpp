#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

// Store expression filenames into a vector 
const vector<string> get_filenames(){
    ifstream file("GEA/hw1/sample_list.txt", ios::in);
    vector<string> filenames;

    if (!file.is_open()){
        std::cout << "\"sample_list.txt\" is failed to be opened!\n";
        // return 1;
    }

    string filename;
    while(getline(file, filename)){
        // cout << filename << endl;
        filenames.push_back(filename);
    }
    file.close();
    return filenames;
}

// Tokenize function for parsing a string into tokens vector
const vector<string> tokenize(string &str, char delimiter){
    vector<string> tokens;
    string token = "";

    // Iterate over characters in the string
    for(int i = 0; i < str.length(); i++){
        if (str[i] != delimiter && i < str.length()-1){
            token += str[i];
        }
        else if (i != str.length()-1){
            tokens.push_back(token);
            token = "";
        }
        else if (str[i] != delimiter && i == str.length()-1){
            tokens.push_back(token+str[i]);
        }
    }
    return tokens;
}


// class BST;
// Define the Nodes
class Node{
    public:
        int entrez;
        vector<int> A_exps; // vector of available gene counts for A cohort
        vector<int> D_exps; // vector of available gene counts for D cohort
        string genename = "";
        float zscore = 0;
        
        // Connect to other nodes
        Node *left;
        Node *right;
        Node *parent;

        // Initialize
        Node():left(0), right(0), parent(0), entrez(0){};
        Node(int ez, string name): left(0), right(0), parent(0), entrez(ez), genename(name){};

        int getez(){return entrez;}
        const vector<int> getAExps(){return A_exps;}
        const vector<int> getDExps(){return D_exps;}

        // friend class BST;
};

class BST{        
    public:
        Node *root = NULL;
        Node* search(int entrez);
        void insert(int entrez, int count, char type, string name);
        void update(Node *current, int count, char type);
        void inorder(Node *current, vector<int>& entrezgenes);
};

Node* BST::search(int entrez){
    Node *current = root;
    while(current != NULL && entrez != current->entrez){
        if (entrez < current->entrez){
            current = current->left;
        } else{
            current = current->right;
        }
    }
    return current;
}

void BST::insert(int entrez, int count, char type, string name){
    Node *x = 0; // Scout
    Node *y = 0; // Pre-parent
    Node *insert_node = new Node(entrez, name); // Name the node
    if (type == 'A'){
        insert_node->A_exps.push_back(count);
    }
    else{
        insert_node->D_exps.push_back(count);
    }

    x = root;
    while (x != NULL){
        y = x;
        if (insert_node->entrez < x->entrez){
            x = x->left;
        } 
        else{
            x = x->right;
        }
    }
    insert_node->parent = y;

    if (y == NULL){
        this->root = insert_node;
    }
    else if (insert_node->entrez < y->entrez){
        y->left = insert_node;
    }
    else{
        y->right = insert_node;
    }
}

void BST::update(Node *current, int count, char type){
    if (type == 'A'){
        current->A_exps.push_back(count);
    }
    else{
        current->D_exps.push_back(count);
    }
}

void BST::inorder(Node *current, vector<int>& entrezgenes){
    if(current == NULL){
        return;
    }
    inorder(current->left, entrezgenes);
    entrezgenes.push_back(current->entrez);
    inorder(current->right, entrezgenes);
}

float mean(vector<int> exps, int n){
    int sum = 0;
    for (int i = 0; i < exps.size(); i++){
        sum += exps[i];
    }
    float mean_value = float(sum)/float(n);
    return mean_value;
}

float variance(vector<int> exps, float mean){
    float var = 0;
    for (int i=0; i < exps.size(); i++){
        float diff = exps[i] - mean;
        var += diff*diff;
    }
    return var;
}

int Partition(vector<float> &v, vector<string> &u, int start, int end){
	int pivot = end;
	int j = start;
	for(int i = start; i < end; ++i){
		if(v[i] < v[pivot]){
			swap(v[i], v[j]);
            swap(u[i], u[j]);
			++j;
		}
	}
	swap(v[j], v[pivot]);
    swap(u[j], u[pivot]);
	return j;
}

void Quicksort(vector<float> &v, vector<string> &u, int start, int end){
	if(start < end){
		int p = Partition(v, u, start, end);
		Quicksort(v, u, start, p-1);
		Quicksort(v, u, p+1, end);
	}
}

int main(void){
    // Build a binary search tree
    BST GeneTree;
    /* Store expression filenames into a vector */
    vector<string> filenames = get_filenames();
    
    // Read gene files iteratively
    // int filecount = 0;
    for(int i=0; i<filenames.size();i++){
        string filename = filenames[i];
        char sampletype = filename[0]; // "A" or "D"

        // Input gene expression file
        ifstream file("GEA/hw1/"+filename, ios::in);
        string geneline;
        if (!file.is_open()){
            std::cout << filenames[i] + " is failed to be opened!\n" << endl;
            return 1;
        }

        // Read lines iteratively
        string line;
        int countline = 0;
        while(getline(file, line)){
            // Skip title line
            if (countline == 0){
                countline++;
                continue;
            }
            // Tokenize the line
            vector<string> tokens = tokenize(line, '\t');
            string gene_name = tokens[1];
            int entrezgene = stoi(tokens[0]);
            int gene_count = stoi(tokens[2]);
            
            // Search if the node has been in tree
            Node *inputGene = GeneTree.search(entrezgene);
            
            if (inputGene == NULL){ // If not found, insert the node
                GeneTree.insert(entrezgene, gene_count, sampletype, gene_name);
            }
            else{ // If found, update the node
                GeneTree.update(inputGene, gene_count, sampletype);
            }
            // Node *requireGene = GeneTree.search(entrezgene);
            countline++;
        }
        file.close();
        // filecount++;
        // cout << "filecount: " << filecount << endl; 
        // if (filecount >= 2){break;}
    }

    // Set initial node to traverse the entire tree inorderly
    Node* initNode = GeneTree.root;
    vector<int> entrezgenes;
    GeneTree.inorder(initNode, entrezgenes);
    std::cout << "Total number of entrez genes: " << entrezgenes.size() << endl;

    // Compute the z-score for each node
    vector<float> z_scores;
    vector<string> gene_names;

    for (int i = 0; i < entrezgenes.size(); i++){
        int entrezgene = entrezgenes[i];
        Node *currentNode = GeneTree.search(entrezgene);
        vector<int> AExps = currentNode->getAExps();
        vector<int> DExps = currentNode->getDExps();
        
        int nA = AExps.size();
        int nD = DExps.size();

        if (nA == 0 || nD == 0){
            continue;
        }

        float meanA = mean(AExps, nA);
        float meanD = mean(DExps, nD);

        float varA = variance(AExps, meanA);
        float varD = variance(DExps, meanD);

        if (varA == 0 || varD == 0){
            continue;
        }

        currentNode->zscore = (meanA-meanD)/sqrt((varA/nA)+(varD/nD));
        z_scores.push_back(currentNode->zscore);
        gene_names.push_back(currentNode->genename);
    }
    
    // Sort the zscore and change orders of gene_names
    Quicksort(z_scores, gene_names, 0, z_scores.size()-1);

    // Write the output file
    ofstream outfile;
    outfile.open("gene_zscores_table_v2.tsv");
    for (int i=0; i<z_scores.size(); i++){
        outfile << gene_names[i] << "\t" << z_scores[i] << endl;
    }
    outfile.close(); 

    return 0;
}
