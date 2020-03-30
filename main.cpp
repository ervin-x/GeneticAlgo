#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int dim, D, PopSize;    // глобальные переменные: dim - количество элементов, D - плотность рюкзака, PopSize - размер попул€ции
double seconds;         // seconds - врем€ работы алгоритма

void ReadCSV()
{
    //from csv
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

void GeneticAlgo(vector<int> Task)
{
    /*
    генераци€ начальной попул€цции
    цикл по i
    TODO VANES
    */
}


void DynamicAlgo(vector<int> Task)
{
    /*
    TODO AZA
    */
}

int main()
{
    ReadCSV();                                              //считываем исходные данные из файла
    for (int k = 1; k <= 100; k++)                          //цикл дл€ генерации задач
    {
        vector<int> Task = TaskGeneration(dim, D);          //Task - текуща€ задача
        clock_t start1 = clock();                           //начало отсчета времени работы алгоритма
        GeneticAlgo(Task);                                  //выполнение задачи методом генетического алгоритма
        clock_t end1 = clock();                             //конец отсчета времени работы алгоритма
        seconds = (double)(end1 - start1) / CLOCKS_PER_SEC; //вычисление времени работы алгоритма
        WriteCSV();                                         //запись результатов работы алгоритма в файл
        clock_t start2 = clock();                           //начало отсчета времени работы алгоритма
        DynamicAlgo(Task);                                  //выполнение задачи методом динамического программировани€
        clock_t end2 = clock();                             //конец отсчета времени работы алгоритма
        seconds = (double)(end2 - start2) / CLOCKS_PER_SEC; //вычисление времени работы алгоритма
        WriteCSV();                                         //запись результатов работы алгоритма в файл
    };
    return 0;
}
