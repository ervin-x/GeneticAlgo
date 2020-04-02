#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>

#define PopSize 10
#define ChromoDim 10
#define NumIterations 100

using namespace std;


void ReadCSV(){
    //from csv
    // TODO MAX
}


void WriteCSV(){
    //to csv
    // TODO MAX
}


vector<int> TaskGeneration(int dim, int D) {
    // возвращает веса элементов vector<int> (dim + int TargetWeight)
    // TODO MAX
}


vector<vector<bool>> PoulationGeneration(int PSize, int dim) {
    /*generate a random population of chromosomes*/

    srand(time(0));
    vector<vector<bool>> Pop(PSize);

    for (auto Сhromo : Pop) {
        for (int i = 0; i < dim; ++i) {
            Сhromo[i] = rand() % 2;
        }
    }

    return Pop;
}


int FitnessFunction(vector<bool> X) {
    // TODO AKELLA
}


vector<vector<bool>> Selection(vector<vector<bool>> Generation) {
    // TODO DIMA
}

vector<bool> SelectionV2(vector<vector<bool>> Generation) {
    // Метод выбора
    // Массив фитнесс функций
    vector<int> chr_fitness;
    for (int i(0); i < Generation.size(); ++i) {
        chr_fitness.push_back(FitnessFunction(Generation[i]));
    }
    // Сортируем массив
    sort(chr_fitness.rbegin(), chr_fitness.rend());
    // Массив индексов хромосом в популяции, отсортированный по убыванию фитнесс функции
    vector<int> indexes;
    for (int i(0); i < chr_fitness.size(); ++i) {
        for (int j(0); j < Generation.size(); ++j) {
            if (chr_fitness[i] == FitnessFunction(Generation[j])) {
                indexes.push_back(j)
            }
        }
    }
    // Выбираем лучшую родительскую хромосому для следующей популяции
    vector<vector<bool>> intermediate_population;
    srand(time(0));
    int groups_count = Generation.size() / 4

    int number = rand() % 100;
    if (number < 50)
        return Generation[indexes[rand() % groups_count]];
    else if ((number >= 50) && (number < 80))
        return Generation[indexes[(rand() % groups_count) + groups_count]];
    else if ((number >= 80) && (number < 95))
        return Generation[indexes[(rand() % groups_count) + 2 * groups_count]];
    else if ((number >= 95) && (number < 100))
        return Generation[indexes[(rand() % groups_count) + 3 * groups_count]];
}

vector<bool> Crossingover(vector<vector<bool>> Generation) {
    // TODO DIMA
}


vector<bool> Mutation(vector<bool> Generation) {
    // TODO AKELLA
}


void GeneticAlgo(vector<int> Task) {
    
    unsigned int start_time =  clock();

    //генерация начальной популяцции
    vector<vector<bool>> InitPop = PoulationGeneration(PopSize, ChromoDim);

    for (int i = 0; i < NumIterations; ++i) {

        if (/*validate on exit*/) {
            //
        }
    }
    
    unsigned int end_time = clock();
    unsigned int work_time = end_time - start_time;
    
}


void DynamicAlgo(vector<int> Task) {
    /*
    TODO AZA
    */
}


int main() {
    // TODO EGOR
}