#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <ctime>
#include <math.h>
#include <time.h>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

#define PopSize 10
#define ChromoDim 10
#define NumIterations 100

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


vector<vector<bool>> PoulationGeneration(int PSize, int dim) {
    /*generate a random population of chromosomes*/

    srand(clock());
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


vector<vector<bool>> Selection(vector<vector<bool>>& Generation)
{
	// TODO DIMA
	// вычисляем сумму фитнесс-функций всех особей
	long long sum_fitness_func(0);
	for (int i(0); i < Generation.size(); ++i)
		sum_fitness_func += FitnessFunction(Generation[i]);

	// вычисляем количество особей i-хромосомы в промежуточной популяции
	vector<int> probabilities_individuals;
	for (int i(0); i < Generation.size(); ++i)
	{
		double hit_probability = FitnessFunction(Generation[i]) * Generation.size() / sum_fitness_func;
		int int_hit_probabilit = round(hit_probability);
		probabilities_individuals.push_back(int_hit_probabilit);
	}

	// формируем промежуточную популяцию
	vector<vector<bool>> intermediate_population;
	for (int i(0); i < probabilities_individuals.size(); ++i)
		for (int j(0); j < probabilities_individuals[i]; ++j)
			intermediate_population.push_back(Generation[i]);

	return intermediate_population;
}


vector<vector<bool>> Crossingover(vector<vector<bool>>& Generation)
{
	// TODO DIMA
	RAND_MAX = Generation.size() - 1;
	srand(time(0));

	vector<vector<bool>> new_population;

	while (new_population.size() != Generation.size())
	{
		int first_parent = rand();
		int second_parent = rand();

		while ((Generation[first_parent] == Generation[second_parent]) &&
			(Generation[first_parent].empty()) &&
			(Generation[second_parent].empty()))
		{
			first_parent = rand();
			second_parent = rand();
		}

		int point_crossingover = rand() % (Generation[first_parent].size() - 2);

		auto child_1 = Generation[first_parent];
		auto child_2 = Generation[second_parent];

		for (int j(point_crossingover); j < Generation[first_parent].size(); ++j)
		{
			child_1[j] = Generation[second_parent][j];
			child_2[j] = Generation[first_parent][j];
		}

		Generation[first_parent].clear();
		Generation[second_parent].clear();

		new_population.push_back(child_1);
		new_population.push_back(child_2);
	}

	return new_population;
}


vector<bool> Mutation(vector<bool> Generation) {
    // TODO AKELLA
}


void GeneticAlgo(vector<int> Task) {
    //TODO VANES
}


vector<int> DynamicAlgo(vector<uint256_t> Task) {
    vector<uint256_t> profits; // стоимости предметов
	vector<uint256_t> weights; // веса предметов
	vector<int> final_set; // результирующий набор предметов
	uint256_t target_weight; // целевой вес
	int subjects = Task.size() - 1; // количество предметов

	// в этой матрице будут храниться решения подзадач данной задачи
	// +1 в размерностях, так как начинаем с 0 как для количества предметов,
	// так и для целевого веса
	uint256_t table[subjects + 1][target_weight + 1];

	for(int i = 0; i < Task.size(); i++) {
		if (i < Task.size() - 1) {
			// векторы введены дополнительно для облегчения решения задачи
			profits.push_back(Task.at(i)); // стоимости равны весам
			weights.push_back(Task.at(i)); // сохраняем веса
			final_set.push_back(0); // изначально ни один предмет не выбран
		} else target_weight = Task.at(i);
	}

	// построение таблицы
	for(int i = 0; i <= subjects; i++) { // для каждого предмета
		for(uint256_t w = 0; w <= target_weight; w++) { // для каждого целевого веса
			// количество предметов и целевой вес = 0
			if (i == 0 || w == 0)
				table[i][w] = 0; // решение задачи очевидно
			else if(weights[i] <= target_weight){
				table[i][w] = max(
					profits[i] + table[i - 1][target_weight - weights[i]],
					table[i - 1][w]
				);
			} else table[i][w] = table[i - 1][w];
		}
	}

	// восстановление отобранных предметов
	uint256_t i = subjects, j = target_weight;
	while (i > 0 && j > 0) {
		if(table[i][j] == table[i - 1][j]){
			// данный предмет не выбран, так как он есть в предыдущей строке
			final_set[i] = 0;
		} else {
			// предмет выбран, так как его не в предыдущей строке
			final_set[i] = 1;
			i--; // переходим к предыдущей строке
			j  = j - weights[i]; // уменьшаем целевой вес
		}
	}

	return final_set;
}


int main() {
    vector<uint256_t> v = TaskGeneration(4, 0.1);
	for (auto i: v){
		cout << i << " ";
	}
	cout << endl;
}