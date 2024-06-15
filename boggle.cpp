#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <random>
#include <chrono>

using namespace std;

#define N 5

unordered_map<string, int> pops;
unordered_map<string, int> states;

int foundstates[50];
int popArr[50];

void loadPops() {
    std::ifstream wordlist("pops.txt");
    string state;
    int population;
    int index = 0;
    while (wordlist >> state >> population) {
        popArr[index++] = population;
        pops[state] = population;
    }
    wordlist.close();

    std::ifstream statelist("states.txt");
    string realstate;
    index = 0;
    while (statelist >> realstate) {
        states[realstate] = index++;
    }
    statelist.close();
    
    for (int i = 0; i < 50; i++) {
        foundstates[i] = 0;
    }
}

void resetFound() {
    for (int i = 0; i < 50; i++) {
        foundstates[i] = 0;
    }
}

struct Trie {
    Trie* children[26];
    bool isEnd;
    Trie() : isEnd(false) {
        for (int i = 0; i < 26; i++)
            children[i] = nullptr;
    }
};

void insert(Trie* root, const string &word) {
    Trie* node = root;
    for (char ch : word) {
        if (!node->children[ch - 'a'])
            node->children[ch - 'a'] = new Trie();
        node = node->children[ch - 'a'];
    }
    node->isEnd = true;
}

bool search(Trie* root, const string &word) {
    Trie* node = root;
    for (char ch : word) {
        if (!node->children[ch - 'a'])
            return false;
        node = node->children[ch - 'a'];
    }
    return node && node->isEnd;
}

bool similarStr(string s1, string s2) {
    int diffs = 0;
    for (int i = 0; i < s1.size(); i++) {
        if (diffs > 1) return false;
        if (s1[i] != s2[i]) diffs++;
    }
    return diffs <= 1;
}

int stateNum(string s) {
    for (auto& entry : states) {
        if (entry.first.size() != s.size()) continue;
        else if (similarStr(entry.first, s)) {
            return entry.second;
        }
    }
    printf("ERROR: couldn't find a state match\n");
    return -1;
}

int getScore(string s) {
    int state_index = stateNum(s);
    if (!foundstates[state_index]) {
        foundstates[state_index] = 1;
        return popArr[state_index];
    }
    else return 0;
}

unordered_set<string> foundWords;
int dx[] = {-1, -1, -1,  0, 0,  1, 1, 1};
int dy[] = {-1,  0,  1, -1, 1, -1, 0, 1};

void solve(int i, int j, vector<vector<char>> &board, Trie* root, string currentWord) {
    if (i < 0 || j < 0 || i >= N || j >= N)
        return;
    char ch = board[i][j];
    if (!root || ch < 'a' || ch > 'z' || !root->children[ch - 'a'])
        return;

    currentWord += ch;
    if (root->children[ch - 'a']->isEnd)
        foundWords.insert(currentWord);
    for (int k = 0; k < 8; k++)
        solve(i + dx[k], j + dy[k], board, root->children[ch - 'a'], currentWord);
}

vector<vector<char>> createBoard(const vector<vector<char>>& cubes, 
                                           const vector<int>& cubePermutations, 
                                           const vector<int>& cubeRotations, 
                                           int size) {
    vector<vector<char>> board(size, vector<char>(size));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int cubeIndex = cubePermutations[i * size + j];
            int faceIndex = cubeRotations[i * size + j];
            board[i][j] = cubes[cubeIndex][faceIndex];
        }
    }

    return board;
}

struct Board {
    vector<vector<char>> cubes;
    vector<int> cubePermutations;
    vector<int> cubeRotations;
    int size;
    int lastScore = 0;
    string unrolled;

    Board(int size, const string& cubeFilename) : size(size) {
        cubes.resize(size * size, vector<char>(6));
        cubePermutations.resize(size * size);
        for (int i = 0; i < size * size; i++) cubePermutations[i] = i;
        cubeRotations.resize(size * size);
        
        random_device rd;
        mt19937 mt(rd());
        uniform_int_distribution<int> dist(0, 5);
        for (int i = 0; i < size * size; i++) cubeRotations[i] = dist(mt);

        ifstream cubefile(cubeFilename);
        if (!cubefile) {
            cerr << "The file " << cubeFilename << " does not exist." << endl;
            return;
        }
        for (int i = 0; i < size * size; i++)
            for (int j = 0; j < 6; j++)
                cubefile >> cubes[i][j], cubes[i][j] += 32;  // Convert to lowercase
    }

    vector<vector<char>> getBoard() {
        vector<vector<char>> board(size, vector<char>(size));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int cubeIndex = cubePermutations[i * size + j];
                int faceIndex = cubeRotations[i * size + j];
                board[i][j] = cubes[cubeIndex][faceIndex];
            }
        }
        return board;
    }
    int evaluateScore(Trie* root, vector<vector<char>> &board) {
        foundWords.clear();
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                solve(i, j, board, root, "");

        vector<string> foundWordsVec(foundWords.begin(), foundWords.end());
        // Sort based on length
        // sort(foundWordsVec.begin(), foundWordsVec.end(), [](const string &a, const string &b) {
        //     return a.length() < b.length();
        // });
        int score = 0;
        for (const auto &w : foundWordsVec)
            // if (w.length() > 3) {
                score += getScore(w);
            // }

        return score;
    }

    void mutate(double temperature, Trie* root) {
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> realDist(0.0, 1.0);
        uniform_int_distribution<int> cubeDist(0, size * size - 1);
        uniform_int_distribution<int> faceDist(0, 5);
        vector<vector<char>> board = getBoard();

        if (realDist(mt) > 0.5) {
            // Switch two random cubes
            vector<int> cubeNums(size * size);
            for (int i = 0; i < size * size; i++) {
                cubeNums[i] = i;
            }
            shuffle(cubeNums.begin(), cubeNums.end(), mt);
            int cube1Pos = cubeNums[0];
            int cube2Pos = cubeNums[1];
            
            // Swap the cubes
            swap(cubePermutations[cube1Pos], cubePermutations[cube2Pos]);
            
            auto board = getBoard();
            int score1 = evaluateScore(root, board); 
            if (realDist(mt) < 1 / (1 + exp(-(score1 - lastScore) / temperature))) {
                lastScore = score1;
            } else {
                // Undo the swap
                swap(cubePermutations[cube1Pos], cubePermutations[cube2Pos]);
            }
        } else {
            // Pick a random cube and rotate it
            int cubeNum = cubeDist(mt);
            int rotation = faceDist(mt);
            int temp = cubeRotations[cubeNum];
            cubeRotations[cubeNum] = rotation;
            auto board = getBoard();
            int score1 = evaluateScore(root, board);
            if (realDist(mt) < 1 / (1 + exp(-(score1 - lastScore) / temperature))) {
                lastScore = score1;
            } else {
                // Undo the rotation
                cubeRotations[cubeNum] = temp;
            }
        }
    }

};

Board runSim(Trie* root, int size, const string& cubeFile) {
    Board board(size, cubeFile);

    double temperature = 1;
    int prevBest = 0;
    int counter = 0;
    while (temperature > 0.0001) {
        board.mutate(temperature * 300, root);  // Pass Trie root to mutate
        if (board.lastScore > prevBest && board.lastScore > 2300) {
            cout << board.lastScore << " ";
            for (const auto& row : board.getBoard()) {
                cout << "{ ";
                for (const auto& ch : row) {
                    cout << ch << " ";
                }
                cout << "} ";
            }
            cout << " ";
            for (const auto& val : board.cubePermutations) cout << val << " ";
            cout << " ";
            for (const auto& val : board.cubeRotations) cout << val << " ";
            cout << temperature << endl;
            prevBest = board.lastScore;
        }
        temperature *= 0.9998;
        counter++;
    }
    return board;
}

int main(int argc, char** argv) {

    Trie* root = new Trie();
    ifstream wordlist("wordlist.txt");
    string word;
    while (wordlist >> word) {
        if(word.length() >= 3)
            insert(root, word);
    }
    wordlist.close();

    const int size = N;  // board size
    const string cubeFile = "cubes.txt";
    int highScore = 0;
    //int foundStates = 0;

    // for (auto& it: pops) {
    //     // Do stuff
    //     cout << it.first << ":\n";
    //     cout << it.second << "\n";
    // }

    loadPops();

    if (argc == 2) {
        vector<vector<char>> board(N, vector<char>(N));

        for (int i = 0; i < N * N; i++) {
            board[i % N][i/N] = (*argv)[i];
        }

        foundWords.clear();
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                solve(i, j, board, root, "");

        vector<string> foundWordsVec(foundWords.begin(), foundWords.end());
        // Sort based on length
        // sort(foundWordsVec.begin(), foundWordsVec.end(), [](const string &a, const string &b) {
        //     return a.length() < b.length();
        // });
        int score = 0;
        for (const auto &w : foundWordsVec)
            // if (w.length() > 3) {
                score += getScore(w);
            // }
        cout << "score: " << score << endl;
        for (int i = 0; i < 50; i++) {
            if (foundstates[i]) cout << "found state #" << i << endl;
        }
        return 0;
    }

    for (int i = 0; i < 10000; i++) {
        resetFound();
        cout << "Starting sim: " << i << ", High score: " << highScore << endl;
        Board b = runSim(root, size, cubeFile);
        cout << "Score: " << b.lastScore << endl;
        cout << b.lastScore << " ";
        for (const auto& row : b.getBoard()) {
            cout << "{ ";
            for (const auto& ch : row) {
                cout << ch << " ";
            }
            cout << "} ";
        }
        cout << " ";
        for (const auto& val : b.cubePermutations) cout << val << " ";
        cout << " ";
        for (const auto& val : b.cubeRotations) cout << val << " ";
        cout << endl;

        if (b.lastScore > highScore) {
            highScore = b.lastScore;
            cout << "New high score: " << highScore << endl;
        }

        for (auto& entry : states) {
            if (foundstates[entry.second]) {
                cout << " found #" << entry.second << ": " << entry.first << endl;
            }
        }
    }

    return 0;   
}
