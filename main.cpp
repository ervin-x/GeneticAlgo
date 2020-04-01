#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

using namespace std;
using namespace chrono;
using namespace boost::multiprecision;

void ReadCSV(){
    //from csv
    // TODO MAX
}


void WriteCSV(){
    //to csv
    // TODO MAX
}


vector<uint256_t> TaskGeneration(int dim, float density) {
    uint256_t max_weight = uint256_t(pow(2, dim / density));
    // cout << "Maxweight = " << max_weight << endl;
    uint256_t sum = 0;
    long long time = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    boost::random::mt19937 mt;
    boost::random::uniform_int_distribution<uint64_t> ui;
    vector<uint256_t> task(dim + 1);
    for (int i = 0; i <= dim; i++) {
        sum += task[i] =  ui(mt) + 1;
        // cout << "Element " << i << " = " << task[i] << endl;;
    }
    // cout << "Sum = " << sum / 2 << endl;
    task[dim] = sum / 2;
    return task;
}

vector<vector<bool>> PopulationGeneration(int PopSize, int dim) {
    // TODO VANES
}


int FitnessFunction(vector<bool> X) {
    // TODO AKELLA
}


vector<vector<bool>> Selection(vector<vector<bool>> Generation) {
    // TODO DIMA
}


vector<bool> Crossingover(vector<vector<bool>> Generation) {
    // TODO DIMA
}


vector<bool> Mutation(vector<bool> Generation) {
    // TODO AKELLA
}


void GeneticAlgo(vector<int> Task) {
    /*
    генерация начальной популяцции
    цикл по i
    TODO VANES
    */
}


void DynamicAlgo(vector<int> Task) {
    /*
    TODO AZA
    */
}


int main() {
    // TODO EGOR
}