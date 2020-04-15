#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <ctime>
#include <sstream>
#include <map>

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

void PrintChromo(const vector<bool>& Chromo) {
	cout << "\nchromo ";
	for (bool gen: Chromo) {
		cout << gen << " ";
	}
}


void PrintGeneration(const vector<vector<bool>>& Generation) {
	for (vector<bool> chromo : Generation) {
		PrintChromo(chromo);
	}
}


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
	cout << "In TaskGeneration" << endl;
    int max_weight = int(pow(2, dim / density));
    //cout << "Maxweight = " << max_weight << endl;
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
	  cout << "In PoulationGeneration" << endl;
    
    srand(clock());
	map <int, bool> int_to_bool = {{1, true}, {0, false}};
    vector<vector<bool>> Pop;

    for (int i = 0; i < PSize; ++i) {
        vector<bool> chromo;
        for (int j = 0; j < dim; ++j) {
            chromo.push_back(int_to_bool[rand() % 2]);
        }
        Pop.push_back(chromo);
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
	if (sum > TWeight) {
		return TWeight;
	} else {
		return TWeight - sum;
	}
}


int FitnessFunctionMod(vector<bool> Chromo, vector<int> Weights, int TWeight) {
	int sum = 0;
	int glob_sum = 0;
	for(int i = 0; i < Chromo.size(); ++i) {
		if (Chromo[i]) {
			sum += Weights[i];
		}
		glob_sum += Weights[i];
	}
	if (sum > TWeight) {
		return 0;
	} else {
		return glob_sum - abs (sum - TWeight);
	}
}


vector<vector<bool>> Selection(const vector<vector<bool>>& Generation, const vector<int>& Weights, const int& TWeight)
{
	cout << "In Selection" << endl;

	srand(clock());

	vector<vector<bool>> intermediate_population;

	while (intermediate_population.size() != Generation.size())
	{
		// количество особей для выбора
		int bidder_counter = 2;
		vector<pair<vector<bool>, int>> bidder_vector;

		int number_bidder = -1;

		for (int i(0); i < bidder_counter; ++i)
		{
			// попытка более-менее разнообразить претендентов на родителей
			while (true)
			{
				int new_number = rand() % (Generation.size());
				if (new_number != number_bidder)
				{
					number_bidder = new_number;
					break;
				}
			}
			cout << "number[" << i << "] = " << number_bidder << "\n";

			int value_fitness_func = FitnessFunction(Generation[number_bidder], Weights, TWeight);
			bidder_vector.push_back(make_pair(Generation[number_bidder], value_fitness_func));
		}

		pair<vector<bool>, int> best_bidder = bidder_vector[0];

		for (int i(1); i < bidder_vector.size(); ++i)
			if (best_bidder.second < bidder_vector[i].second)
				best_bidder = bidder_vector[i];

		cout << "Best biddet " << best_bidder.second << "\n";

		intermediate_population.push_back(best_bidder.first);
	}

	return intermediate_population;
}

vector<vector<bool>> Crossingover(vector<vector<bool>>& Generation)
{
	cout << "In Crossingover" << endl;
	srand(clock());

	vector<vector<bool>> new_population;

	while (new_population.size() != Generation.size())
	{
		cout << "In Crossingover 1" << endl;

		int first_parent = -1;
		int second_parent = -1;

		while (true)
		{
			first_parent = rand() % (Generation.size());
			second_parent = rand() % (Generation.size());

			cout << "first_parent " << first_parent << "second_parent " << second_parent << endl;

			if (Generation[first_parent].size() != 0 &&
				Generation[second_parent].size() != 0 &&
				first_parent != second_parent)
				break;
			else
				continue;
		}

		int point_crossingover = rand() % (Generation[first_parent].size() - 1);

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
	cout << "In Mutation" << endl;
	// настраиваем ГПСЧ
	long long time = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
    	boost::random::mt19937 mt(time);
    	boost::random::uniform_real_distribution<double> ui(0, 1);
	// пробегаемся по всему поколению
	for(int i = 0; i < Generation.size(); i++){
	        // случайным образом выбираем ген для инвертирования
	        int genToChange = rand() % Generation[i].size();
	        // с вероятностью MutationChance инвертируем данный ген
	        if (MutationChance > ui(mt))
			Generation[i][genToChange] = !Generation[i][genToChange];
	}
		
	return Generation;
}


vector<bool> GeneticAlgo(const vector<int>& Task, const int& PSize, const int& NumIterations) {

	cout << "\nIn GeneticAlgo" << endl;
	int dim = Task.size() - 1;
	int TWeight = Task.back();
	vector<int> Weights(Task.begin(), Task.end()-1);
    vector<vector<bool>> population = PoulationGeneration(PSize, dim);

	//cout << "population.size() = " << population.size() << endl;

	for (int i = 0; i < NumIterations; ++i) {

		population = Selection(population, Weights, TWeight);
		population = Crossingover(population);
		population = Mutation(population);

		int fit_value = TWeight;
		int repeat_num = 0;
		vector<bool> best_chromo;

		// выбор лучшей хромосомы по фитнес функции
		for (auto chromo : population) {
			int _fit_value = FitnessFunctionMod(chromo, Weights, TWeight);
			if (_fit_value > fit_value) {
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

	// в этой матрице будут храниться решения подзадач данной задачи
	// +1 в размерностях, так как начинаем с 0 как для количества предметов,
	// так и для целевого веса
	int** table = new int*[subjects + 1];
	for(int i = 0; i < subjects + 1; i++)
		table[i] = new int[target_weight + 1];

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

	// очистка памяти
	for(int i = 0; i < subjects + 1; i++)
		delete table[i];
	delete[] table;

	return final_set;
}

int main()
{
	int num_file_lines = 4;

	// int line_number = 0;
	float density;
	int pop_size;
	float mutation_prob;

	for (int l = 0; l < num_file_lines; ++l) {

		ReadCSV(ifilename, l, density, pop_size, mutation_prob);
		cout << "input parameters:"
			<< " line_number: " << l
			<< " density: " << density
			<< " pop_size: " << pop_size
			<< " mutation_prob: " << mutation_prob
			<< endl;

		long long gen_time = 0;
		long long dynamic_time = 0;

		for (int k = 1; k <= 100; k++)
		{
			cout << "k = " << k << endl;

			vector<int> Task = TaskGeneration(pop_size, density);

			clock_t start_time = clock();
			vector<bool> bestGAChromo = GeneticAlgo(Task, pop_size, NumIterations);
			gen_time += (long long)(clock() - start_time) / CLOCKS_PER_SEC;

			start_time = clock();
			vector<int> bestDynChromo = DynamicAlgo(Task);
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

	}

    return 0;
}
