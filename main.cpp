#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <ctime>
#include <math.h>
#include <time.h>
#include <sstream>

#include <boost/random.hpp>

using namespace std;
using namespace chrono;

const string ifilename = "in.csv";
const string ofilename = "out.csv";

const int PopSize = 10;
const int ChromoDim = 10;
const int NumIterations = 100;
const int ExitThreshold = 10;
const float MutationChance = 0.001; // вероятность мутации


void ReadCSV(const string& ifilename, int line_number, float &density, int &pop_size, float &mutation_prob) {
    ifstream file(ifilename);
    if(!file.is_open()) throw runtime_error("Could not open input file");
    string line, substr;
    for (int i = 0; i < line_number + 1; i++) {
        getline(file, line);
    }
    stringstream ss(line);
    
    getline(ss, substr, ',');
    density = stof(substr);

    getline(ss, substr, ',');
    pop_size = stoi(substr);

    getline(ss, substr, ',');
    mutation_prob = stof(substr);
    file.close();
}

  
void WriteCSV(const string& ofilename, float density, int dim, long long gen_time, long long dynamic_time){
    ofstream file(ofilename, fstream::app);
    file << density << ',' << dim << ',' << gen_time << ',' << dynamic_time << endl;
    file.close();
}


vector<int> TaskGeneration(int dim, float density) {
    int max_weight = int(pow(2, dim / density));
    // cout << "Maxweight = " << max_weight << endl;
    int sum = 0;
    long long time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
    boost::random::mt19937 mt(time);
    boost::random::uniform_int_distribution<int> ui(int(1), max_weight);
    vector<int> task(dim + 1);
    for (int i = 0; i <= dim; i++) {
        sum += task[i] =  ui(mt) + 1;
        // cout << "Element " << i << " = " << task[i] << endl;;
    }
    // cout << "Sum = " << sum / 2 << endl;
    task[dim] = sum / 2;
    return task;
}


vector<vector<bool>> PoulationGeneration(int PSize, int dim) {

    srand(clock());
    vector<vector<bool>> Pop(PSize);

    for (auto Chromo : Pop) {
        for (int i = 0; i < dim; ++i) {
            Chromo.push_back(rand() % 2);
        }
    }

    return Pop;
}


int FitnessFunction(vector<bool> Chromo, vector<int> Weights, int TWeight) { // DONE -> AZAMAT
	int sum = 0;
	for(int i = 0; i < Chromo.size(); ++i) {
		if (Chromo[i]) {
			sum += Weights[i];
		}
	}
	return abs (TWeight - sum);
}


vector<vector<bool>> Selection(vector<vector<bool>>& Generation, vector<int> Weights, int TWeight)
{
	// вычисляем сумму фитнесс-функций всех особей
	long long sum_fitness_func(0);
	for (int i(0); i < Generation.size(); ++i)
		sum_fitness_func += FitnessFunction(Generation[i], Weights, TWeight);

	// вычисляем количество особей i-хромосомы в промежуточной популяции
	vector<int> probabilities_individuals;
	for (int i(0); i < Generation.size(); ++i)
	{
		double hit_probability = FitnessFunction(Generation[i], Weights, TWeight) * Generation.size() / sum_fitness_func;
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
	// RAND_MAX = Generation.size() - 1;
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


vector<vector<bool>> Mutation(vector<vector<bool>> Generation) { // DONE -> AZAMAT
	// настраиваем ГПСЧ
	long long time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
    	boost::random::mt19937 mt(time);
    	boost::random::uniform_real_distribution<double> ui(0, 1);
	// случайным образом выбираем ген для инвертирования
	int chronoToChange = rand() % Generation.size();
	// с вероятностью MutationChance инвертируем данный ген
	if (MutationChance > ui(mt))
		for (bool Gen: Generation[chronoToChange]){
			Gen = !Gen;
		}

	return Generation;
}


vector<bool> GeneticAlgo(const vector<int>& Task, const int& PSize, const int& NumIterations) {

	int dim = Task.size() - 1;
	int TWeight = Task.back();
	vector<int> Weights(Task.begin(), Task.end()-1);
    vector<vector<bool>> population = PoulationGeneration(PSize, dim);

	for (int i = 0; i < NumIterations; ++i) {

		population = Selection(population, Weights, TWeight);
		population = Crossingover(population);
		population = Mutation(population);

		int fit_value = TWeight;
		int repeat_num = 0;
		vector<bool> best_chromo;

		// выбор лучшей хромосомы по фитнес функции
		for (auto chromo : population) {
			int _fit_value = FitnessFunction(chromo, Weights, TWeight);
			if (_fit_value < fit_value) {
				fit_value = _fit_value;
				best_chromo = chromo;
			} else if (_fit_value == fit_value) {
				if (repeat_num >= ExitThreshold) {

				return best_chromo;

				} else {
					++repeat_num;
				}
			} 
		}
	}
}

vector<int> DynamicAlgo(vector<int> Task) { // DONE -> AZAMAT
	vector<int> profits; // стоимости предметов
	vector<int> weights; // веса предметов
	vector<int> final_set; // результирующий набор предметов
	int target_weight; // целевой вес
	int subjects = Task.size() - 1; // количество предметов

	// в этой матрице будут храниться решения подзадач данной задачи
	// +1 в размерностях, так как начинаем с 0 как для количества предметов,
	// так и для целевого веса
	int table[subjects + 1][target_weight + 1];

	// предмет с нулевым весом и стоимостью
	// необходим для общности решения
	profits.push_back(0);
	weights.push_back(0);

	for(int i = 0; i < Task.size(); i++) {
		if (i < Task.size() - 1) {
			// векторы введены дополнительно для облегчения решения задачи
			profits.push_back(Task.at(i)); // стоимости равны весам
			weights.push_back(Task.at(i)); // сохраняем веса
			final_set.push_back(0); // изначально ни один предмет не выбран
		}
		else target_weight = Task.at(i); // сохраняется целевой вес
	}

	// построение таблицы
	for(int i = 0; i <= subjects; i++) { // для каждого предмета
		for(int w = 0; w <= target_weight; w++) { // для каждого целевого веса
			// количество предметов и целевой вес = 0
			if (i == 0 || w == 0)
				table[i][w] = 0; // решение задачи очевидно
			else if(weights[i] <= w){ // вес предмета меньше или равен рассматриваемому целевому весу
				table[i][w] = max(
					profits[i] + table[i - 1][w - weights[i]],
					table[i - 1][w]
				);
			} else table[i][w] = table[i - 1][w]; // вес предмета больше рассматриваемого целевого веса
		}
	}

	// восстановление отобранных предметов
	int i = subjects;
	int j = target_weight;
	while (i > 0 && j > 0) {
		if(table[i][j] == table[i - 1][j]){
			// данный предмет не выбран, так как он есть в предыдущей строке
			final_set[i - 1] = 0;
		} else {
			// предмет выбран, так как его нет в предыдущей строке
			final_set[i - 1] = 1;
			j = j - weights[i]; // уменьшаем целевой вес для продолжения поиска
		}
		i--; // переходим к предыдущей строке
	}

	return final_set;
}

int main()
{
	int line_number = 0;
	float density;
	int pop_size;
	float mutation_prob;

    ReadCSV(ifilename, line_number, density, pop_size, mutation_prob);
	cout << "input parameters:"
		 << " line_number: " << line_number
		 << " density: " << density
		 << " pop_size: " << pop_size
		 << " mutation_prob: " << mutation_prob
		 << endl;

	long long gen_time = 0;
	long long dynamic_time = 0;

    for (int k = 1; k <= 100; k++)
    {
        vector<int> Task = TaskGeneration(pop_size, density);

        clock_t start_time = clock();
        GeneticAlgo(Task, pop_size, NumIterations);
        gen_time += (long long)(clock() - start_time) / CLOCKS_PER_SEC;

        start_time = clock();
        DynamicAlgo(Task);
        dynamic_time += (long long)(start_time - start_time) / CLOCKS_PER_SEC;

    }

	gen_time /= 100;
	dynamic_time /= 100;

	WriteCSV(ofilename, density, pop_size, gen_time, dynamic_time);
	cout << "output parameters:"
		 << " density: " << density
		 << " pop_size: " << pop_size
		 << " gen_time: " << gen_time
		 << " dynamic_time: " << dynamic_time
		 << endl;

    return 0;
}
