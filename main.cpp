#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#define PopSize 10
#define ChromoDim 10
#define NumIterations 100

using namespace std;


void ReadCSV() {
    // TODO MAX
}


void WriteCSV()
{
    //to csv
    // TODO MAX
}

vector<int> TaskGeneration(int dimension, int Density)
{
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


vector<bool> Crossingover(vector<vector<bool>> Generation) {
    // TODO DIMA
}


vector<bool> Mutation(vector<bool> Generation) {
    // TODO AKELLA
}


void GeneticAlgo(vector<int> Task) {

    //генерация начальной популяцции
    vector<vector<bool>> InitPop = PoulationGeneration(PopSize, ChromoDim);

    for (int i = 0; i < NumIterations; ++i) {

        if (/*validate on exit*/) {
            //
        }
    }
}


int main()
{
    ReadCSV();                                              //считываем исходные данные из файла
    for (int k = 1; k <= 100; k++) {                          //цикл для генерации задач
        vector<int> Task = TaskGeneration(dim, D);          //Task - текущая задача
        clock_t start1 = clock();                           //начало отсчета времени работы алгоритма
        GeneticAlgo(Task);                                  //выполнение задачи методом генетического алгоритма
        clock_t end1 = clock();                             //конец отсчета времени работы алгоритма
        seconds = (double)(end1 - start1) / CLOCKS_PER_SEC; //вычисление времени работы алгоритма
        WriteCSV();                                         //запись результатов работы алгоритма в файл
        clock_t start2 = clock();                           //начало отсчета времени работы алгоритма
        DynamicAlgo(Task);                                  //выполнение задачи методом динамического программирования
        clock_t end2 = clock();                             //конец отсчета времени работы алгоритма
        seconds = (double)(end2 - start2) / CLOCKS_PER_SEC; //вычисление времени работы алгоритма
        WriteCSV();                                         //запись результатов работы алгоритма в файл
    };
    return 0;

}
