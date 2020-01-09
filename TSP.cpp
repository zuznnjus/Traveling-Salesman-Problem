#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

#define ALPHA 9
#define BETA 12
#define Q_0 0.2
#define EVAPORATION 0.4
#define ITERATIONS 200
#define DIVISOR 50
#define NUMBER_OF_ANTS 50
#define INITIAL_PHEROMONE 0.001

struct Point {
	double x;
	double y;
	bool visited;
};

struct Edge {
	double distance;
	double pheromone;
	double eta;
};

void generateFile() {
	srand(time(NULL));
	ofstream file("plik100.txt");
	int size = 100,x,y;
	file << size<<endl;
	for (int i = 0; i < size; i++) {
		x = rand() % 1000 + 1;
		y = rand() % 1000 + 1;
		file << i + 1 << " " << x << " " << y << endl;
	}
	file.close();
}

Point* readFile(const string &fileName, int& size)
{
	ifstream file(fileName);
	if (file.is_open())
		file >> size;

	int tmp;
	Point* coordinates = new Point[size];

	for (int i = 0; i < size; i++)
	{
		file >> tmp;
		file >> coordinates[i].x;
		file >> coordinates[i].y;
		coordinates[i].visited = false;
	}
	return coordinates;
}

double calculateDistance(Point point1, Point point2)
{
	double distance = sqrt((point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y));
	return distance;
}

Edge** distances(Point* coordinates, const int &size) {
	Edge** edges = new Edge * [size];
	for (int i = 0; i < size; ++i)
		edges[i] = new Edge[size];

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i < j) {
				edges[i][j].distance = calculateDistance(coordinates[i], coordinates[j]);
				edges[i][j].eta = 1 / edges[i][j].distance;
				edges[i][j].pheromone = INITIAL_PHEROMONE;
			}
		}
	}
	return edges;
}

double greedy(Point* coordinates, Edge** edges, const int &size) {
	double* distances = new double[size];
	int a, b;

	for (int point = 1; point < size + 1; point++) {
		double min = 0, distance, sum = 0;
		int pointNumber = point, nextPointNumber = 0, index = 1;

		coordinates[pointNumber - 1].visited = true;
		while (index < size) {

			for (int i = 0; i < size; i++) {
				if (!coordinates[i].visited) {
					((pointNumber - 1) < i) ? (a = pointNumber - 1, b = i) : (a = i, b = pointNumber - 1);
					min = edges[a][b].distance;
					nextPointNumber = i + 1;
					break;
				}
			}

			for (int i = 0; i < size; i++) {
				if (!coordinates[i].visited) {
					((pointNumber - 1) < i) ? (a = pointNumber - 1, b = i) : (a = i, b = pointNumber - 1);
					distance = edges[a][b].distance;
					if (distance < min) {
						min = distance;
						nextPointNumber = i + 1;
					}
				}
			}
			coordinates[nextPointNumber - 1].visited = true;
			sum += min;
			pointNumber = nextPointNumber;
			index++;
		}

		((pointNumber - 1) < (point - 1)) ? (a = pointNumber - 1, b = point - 1) : (a = point - 1, b = pointNumber - 1);
		sum += edges[a][b].distance;
		distances[point - 1] = sum;

		for (int i = 0; i < size; i++)
		{
			coordinates[i].visited = false;
		}
	}
	double minDist = distances[0];
	for (int i = 1; i < size; i++)
	{
		if (distances[i] < minDist) {
			minDist = distances[i];
		}
	}
	return minDist;
}

int probability(Edge** edges, const int &currentCity, bool* visitedCities, const int &size) {
	double sum = 0, prob, probMax;
	double* numerators = new double[size];
	int nextCity = 0, a, b;

	for (int i = 0; i < size; i++)
	{
		if (!visitedCities[i]) {
			(currentCity < i) ? (a = i, b = currentCity) : (a = currentCity, b = i);
			sum += pow(edges[b][a].pheromone, ALPHA) * pow(edges[b][a].eta, BETA);
			numerators[i] = pow(edges[b][a].pheromone, ALPHA) * pow(edges[b][a].eta, BETA);
		}
		else
			numerators[i] = -1;
	}

	for (int i = 0; i < size; i++) {
		if (numerators[i] > -1) {
			nextCity = i;
			probMax = numerators[i] / sum;
		}
	}

	for (int i = 0; i < size; i++) {
		if (numerators[i] > -1) {
			prob = numerators[i] / sum;
			if (prob > probMax) {
				probMax = prob;
				nextCity = i;
			}
		}
	}
	delete[]numerators;
	return nextCity;
}

void updatePheromoneOnlyBestAnt(Edge** edges, int** antsRoutes, double* antsDistances, const int &size,const int &bestAnt) {
	double addedValue = Q_0 / antsDistances[bestAnt];
	int a, b;
	for (int i = 0; i < size; i++) {
		(antsRoutes[bestAnt][i] < antsRoutes[bestAnt][i + 1]) ? (a = antsRoutes[bestAnt][i], b = antsRoutes[bestAnt][i + 1]) : (a = antsRoutes[bestAnt][i + 1], b = antsRoutes[bestAnt][i]);
		edges[a][b].pheromone += addedValue;
	}
}

void evaporatePheromone(Edge** edges, const int &size) {
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i < j)
				edges[i][j].pheromone = (1 - EVAPORATION) * edges[i][j].pheromone;
		}
	}
}

double antColony(Edge** edges, Point* cooridnates, const int &size) {
	int nextPoint, bestAnt, repetitionCounter = 0, a, b;
	double localBestResult, bestResult;
	double* arrayLocalBestResult = new double[ITERATIONS];
	srand(time(NULL));

	for (int it = 0; it < ITERATIONS; it++) {

		int** antsRoutes = new int* [NUMBER_OF_ANTS];
		for (int i = 0; i < NUMBER_OF_ANTS; ++i)
			antsRoutes[i] = new int[size + 1];

		double* antsDistances = new double[NUMBER_OF_ANTS];
		for (int i = 0; i < NUMBER_OF_ANTS; ++i)
			antsDistances[i] = 0;

		for (int ant = 0; ant < NUMBER_OF_ANTS; ant++) {

			bool* visitedCities = new bool[size];
			for (int i = 0; i < size; i++) {
				visitedCities[i] = false;
			}

			int startPoint = rand() % size;
			int currentPoint = startPoint;
			visitedCities[currentPoint] = true;
			antsRoutes[ant][0] = currentPoint;

			int index = 1;
			while (index < size) {
				nextPoint = probability(edges, currentPoint, visitedCities, size);
				antsRoutes[ant][index] = nextPoint;
				visitedCities[nextPoint] = true;
				currentPoint = nextPoint;
				index++;
			}
			antsRoutes[ant][index] = startPoint;
		}

		for (int ant = 0; ant < NUMBER_OF_ANTS; ant++)
		{
			for (int i = 0; i < size; i++) {
				(antsRoutes[ant][i] < antsRoutes[ant][i + 1]) ? (a = antsRoutes[ant][i], b = antsRoutes[ant][i + 1]) : (a = antsRoutes[ant][i + 1], b = antsRoutes[ant][i]);
				antsDistances[ant] += edges[a][b].distance;
			}
			if (ant == 0) {
				localBestResult = antsDistances[ant];
				bestAnt = ant;
			}
			else
				if (localBestResult > antsDistances[ant]) {
					localBestResult = antsDistances[ant];
					bestAnt = ant;
				}
		}
		cout << localBestResult << endl;

		arrayLocalBestResult[it] = localBestResult;
		if (it > 0) {
			if (int(arrayLocalBestResult[it - 1]) == int(arrayLocalBestResult[it]))
				repetitionCounter++;
			else
				repetitionCounter = 0;
		}

		if (repetitionCounter > (ITERATIONS / DIVISOR)) {
			cout << "---wygladzenie krawedzi---" << endl;
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					if (i < j) {
						edges[i][j].pheromone = INITIAL_PHEROMONE * (sqrt(edges[i][j].pheromone / INITIAL_PHEROMONE));
					}
				}
			}
		}

		if (it == 0) {
			bestResult = localBestResult;
		}
		if (bestResult > localBestResult) {
			bestResult = localBestResult;
		}

		evaporatePheromone(edges, size);
		updatePheromoneOnlyBestAnt(edges, antsRoutes, antsDistances, size, bestAnt);

		delete[]antsDistances;
		delete[]antsRoutes;
	}

	delete[]arrayLocalBestResult;
	return bestResult;
}

int main()
{
	int size = 1;
	string fileName = "plik100.txt";
	fstream file;
	file.open("wyniki.txt", ios::out | ios::app);
	Point* coordinates = readFile(fileName, size);
	Edge** edges = distances(coordinates, size);
	double bestResultGreedy = greedy(coordinates, edges, size);
	clock_t begin = clock();
	double bestResultACO = antColony(edges, coordinates, size);
	clock_t end = clock();
	cout << "best greedy result: " << bestResultGreedy << endl;
	cout << "best ACO result: " << bestResultACO << endl;
	cout << "ACO time: " << (float)(end - begin) / CLOCKS_PER_SEC<<" s"<<endl;
	file << bestResultGreedy << ";  " << bestResultACO << endl;
	file.close();
	delete[]edges;
	delete[]coordinates;
	//generateFile();
	return 0;
}