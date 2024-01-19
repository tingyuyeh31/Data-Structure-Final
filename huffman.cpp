#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

using namespace std;

const map<char, int> get_freq(){
    // Read file
    ifstream file("pdb_seqres/pdb_seqres.txt", ios::in);
    if (!file.is_open()){
        cout << "Error: \"pdb_seqres.txt\" not found" << endl;
    }

    // Create a hash map ("A":frequency)
    map<char, int> freq;
    string line;
    int countline = 0;
    while(getline(file, line)){
        if (line[0] == '>'){ // Skip identifier lines
            countline++;
            continue;
        }
        else{
            for (int i = 0; i < line.length(); i++){
                char letter = line[i];
                if (!freq.count(letter)){
                    freq[letter] = 1;
                }
                else{
                    freq[letter]++;
                }
                // cout << letter << " " << freq.count(letter) << " " << freq[letter] << endl;
            }
            countline++;
        }
    }
    return freq;
}

// define object "node"
class node{
	public:
		node(){amount=0; letter='\0'; left=NULL; right=NULL;};
		char letter;
		int amount;
		node *left;
		node *right;
};

// Sorted by amounts
void bubblesort(vector<node> &huffman){
    for (int j = 0; j < huffman.size(); j++){
        for (int i = 0; i < huffman.size()-1; i++){
            if (huffman[i].amount > huffman[i+1].amount){
                node temp = huffman[i];
                huffman[i] = huffman[i+1];
                huffman[i+1] = temp;
            }
        }
    }
}

// Construct Huffman Tree
// Like heap structure
void huffman_tree(vector<node> &huffman){
    while(huffman.size() != 1){
        node temp;
        temp.amount = huffman[0].amount + huffman[1].amount;
        node *L = new node;
        node *R = new node;

        L->letter = huffman[0].letter;
        L->amount = huffman[0].amount;
        L->left=huffman[0].left;
        L->right=huffman[0].right;

        R->letter = huffman[1].letter;
        R->amount = huffman[1].amount;
        R->left=huffman[1].left;
        R->right=huffman[1].right;
        temp.left = L;
        temp.right = R;

        huffman.erase(huffman.begin(),huffman.begin()+2); // Remove index 0 and 1 in huffman
        
        // Find the location to insert the merged node (temp)
        int i = 0;
        while(huffman[i].amount < temp.amount && i < huffman.size()){
            i++;
        }
        huffman.insert(huffman.begin()+i, temp);
    }
}

// recode the huffman code for each letter (AA)
void recursive(node *huffman, map<char, string> &code_map, string &code){
    if (huffman->left != NULL && huffman->right != NULL){
        code += "0";
        recursive(huffman->left, code_map, code);
        code += "1";
        recursive(huffman->right, code_map, code);
        code.pop_back();
    }
    else{
        code_map[huffman->letter] = code;
        code.pop_back();
        return;
    }
}

void encode(map<char, string> &code_map){
    ifstream infile("pdb_seqres/pdb_seqres.txt", ios::in);
    ofstream outfile;
    outfile.open("huffman_encoded.txt");
    int countline = 0;
    string line;
    while(getline(infile, line)){
        if (line[0] == '>'){
            countline++;
            outfile << line << endl;
            continue;
        }
        else{
            for (int i = 0; i < line.length(); i++){
                char letter = line[i];
                string bitstring = code_map[letter];
                outfile << bitstring;
            }
            outfile << endl;
            countline++;
        }
    }
    outfile.close();
}

void decode(node *huffman){
    ifstream infile("huffman_encoded.txt", ios::in);
    ofstream outfile;
    outfile.open("huffman_decoded.txt");
    int countline = 0;
    string line;
    node *root = huffman;
    while(getline(infile, line)){
        if (line[0] == '>'){
            countline++;
            outfile << line << endl;
            continue;
        }
        else{
            string decode_string;
            for (int i = 0; i < line.length(); i++){
                char bit = line[i];
                if (bit == '0'){
                    huffman = huffman->left;
                }else{
                    huffman = huffman->right;
                }
                if (huffman->letter != '\0'){
                    char decodeAA = huffman->letter;
                    decode_string.push_back(decodeAA);
                    huffman = root;
                }
            }
            countline++;
            outfile << decode_string << endl;
        }
    }
    outfile.close();
}

int main(void){
    cout << "Frequency counting!" << endl;
    map<char, int> freq = get_freq();
    map<char, int>::iterator iter;

    // Save nodes into a vector
    cout << "Leaf nodes creating!" << endl;
    vector<node> huffman;
    // vector<int> AAnames;
    for (iter=freq.begin(); iter!=freq.end(); iter++){
        node temp;
        temp.letter = iter->first;
        temp.amount = iter->second;
        huffman.push_back(temp);
        // AAnames.push_back(iter->first);
        // cout << iter->first << ": " << iter->second << endl;
    }

    // Bubble sorting with amounts
    cout << "Bubble sorting!" << endl;
    bubblesort(huffman);

    // Construct Huffman Tree
    cout << "Huffman tree constructing!" << endl;
    huffman_tree(huffman);

    // Find the leaf nodes that contain the letters and save their own codes
    cout << "Bit code table creating!" << endl;
    map<char, string> code_map; // {"A": 1110, ...}
    string code = "";
    recursive(&huffman[0], code_map, code);

    // Print the coding table
    map<char, string>::iterator iter2;
    ofstream outfile;
    outfile.open("AA_bitstring_table.txt");
    for (iter2=code_map.begin(); iter2!=code_map.end(); iter2++){
        outfile << iter2->first << "\t" << iter2->second << endl;
    }
    outfile.close();

    cout << "Encoding!" << endl;
    encode(code_map);
    cout << "Decoding!" << endl;
    decode(&huffman[0]);

    cout << "Done!" << endl;
    return 0;
}