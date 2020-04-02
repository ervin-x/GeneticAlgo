#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <math.h>
#include <time.h>

#define PopSize 10
#define ChromoDim 10
#define NumIterations 100
#define ExitThreshold 10

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


vector<bool> Crossingover(vector<vector<bool>>& Generation)
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


vector<bool> GeneticAlgo(vector<int> Task, int PSize, int dim, int NumIterations) {

    vector<vector<bool>> population = PoulationGeneration(PSize, dim);

	for (int i = 0; i < NumIterations; ++i) {

		population = Selection(population);
		population = Crossingover(population);
		population = Mutation(population);

		int fit_value = 0;
		int repeat_num = 0;
		vector<bool> best_chromo;

		// выбор лучшей хромосомы по фитнес функции
		for (auto chromo : population) {

			int _fit_value = FitnessFunction(chromo);

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


int main() {
    // TODO EGOR
}