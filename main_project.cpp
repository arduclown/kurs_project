#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
using namespace std;

struct Rectangle
{
    int position;
    int x_left;
    int x_right;
    int y_down;
    int y_up;
};

double e = 1E-13; //невязка
int m_ieration = 1000; //максимальное количество итераций
int cnt_iter = 0; //счетчик итераций
int x_size = 0; //кол-во элементов по оси х
int y_size = 0; //кол-во элементов по оси у
int global_size = 0; //общее кол-во элементов
double h_x = 0; //шаг сетки по х
double h_y = 0; //шаг сетки по у
int L = 0; //кол-во подобластей
int cntKraev = 0; //кол-во краевых условий

vector<int> X_grid; //разбиение на сетку по х (координаты)
vector<int> Y_grid; //разбиение на сетку по у (координаты)
vector<vector<double> > G; //матрица жесткости
vector<vector<double> > M; //матрица масс
vector<double> b; //вектор правых частей
vector<Rectangle> rectangles;
vector<vector<int> > kraevye_uslov; //кравеые условия

//хранение матрицы в разреженном строчном формате
vector<double> di; //диагональ
vector<double> gg; //ненулевые компоненты матрицы
vector<int> ig; //индексы строк ненулевых компонет
vector<int> jg; //индексы столбцов ненулевых компонент

double lambdaV (double x, double y)
{
    return 1;
}

double gammaV(double x, double y)
{
    return 1;
}

int NumberOfUnknowns(int elem) { return 4; }

int IndexOfUnknown(int ielem, int localIndex, int Nx) {
    // Примерный способ определения глобального индекса узла
    // ielem - индекс конечного элемента
    // localIndex - локальный индекс узла внутри элемента (0, 1, 2 или 3)
    // Nx - количество узлов по оси x

    // Предположим, что элементы нумеруются последовательно по строкам (слева направо, сверху вниз)
    int elementsPerRow = Nx - 1; // Количество элементов в строке
    int row = ielem / elementsPerRow; // Определение строки
    int col = ielem % elementsPerRow; // Определение столбца

    // Возвращаем глобальный индекс узла на основе его локального индекса
    switch (localIndex) {
        case 0: // Левый нижний узел
            return row * Nx + col;
        case 1: // Правый нижний узел
            return row * Nx + col + 1;
        case 2: // Левый верхний узел
            return (row + 1) * Nx + col;
        case 3: // Правый верхний узел
            return (row + 1) * Nx + col + 1;
        default:
            cerr << "Invalid local index for element" << endl;
            return -1; // Ошибка, если индекс выходит за пределы допустимого диапазона
    }
}




//построение сетки
void GridBuilder()
{
    ifstream input("matrix.txt");
    input >> x_size >> y_size; 
    X_grid.resize(x_size);
    Y_grid.resize(y_size);
    for (int i = 0; i < x_size; i++)
        input >> X_grid[i];
    for (int i = 0; i < y_size; i++)
        input >> X_grid[i];
    input.close();
    global_size = x_size * y_size;

    //построение подобластей    
    input >> L;
    rectangles.resize(L);
    for (int i = 0; i < L; i++)
        input >> rectangles[i].position >> rectangles[i].x_left >> rectangles[i].x_right >> rectangles[i].y_down >> rectangles[i].y_up;
    input.close();
}





//построение краевых условий
void Kraevye() //описание: 1 - тип кравевого усл, 2 - номер границы, 3 - координата x начало, 4 - координата x конец, 5 - координата y начало, 6 - координата y конец
{
    ifstream input("kraevie.txt");
    input >> cntKraev;
    kraevye_uslov.resize(cntKraev);
    for (int i = 0; i < cntKraev; i++)
    {
        kraevye_uslov[i].resize(6);
        for (int j = 0; j < 6; j++)
        {
            input >> kraevye_uslov[i][j];
        }
    }
    input.close();
}

void Third_K(vector<int> th_krU)
{
    //параллельно у
    if (th_krU[2] == th_krU[3])
    {

    }
    else
    {

    }
}

void Second_K(vector<int> se_krU)
{

}

void First_K(vector<int> fi_krU)
{

}



//портрет матрицы. пока не работает! доделать!!!!!!
void Portret_builders() //построение портрета
{
    ig.resize(global_size + 1);
    vector<set<int> > list; // временный массив для хранения списка
    list.resize(global_size);

    for (int ielem = 0; ielem < L; ielem++) //идем по всем конечным элементам
    {
        vector<int> nodes;
        for (int i = 0; i < NumberOfUnknowns(ielem); ++i) {
            nodes.push_back(IndexOfUnknown(ielem, i, x_size));
        }
        
        //формируем связь между элементами
        for (int i = 0; i < NumberOfUnknowns(ielem); i++)
        {
            
            for (int j = i + 1; j < NumberOfUnknowns(ielem); j++)
            {
                int ind1 = min(nodes[i], nodes[j]);
                int ind2 = max(nodes[i], nodes[j]);
                list[ind2].insert(ind1);

            }
        }
    }

    ig[0] = 0;
    for (int i = 0; i < global_size; i++)
    {
        for (int index: list[i]) jg.push_back(index);
        ig[i+1] = jg.size();
    }

}

void LocalG_matrix(double x, double y) //построение локальной матрицы жесткости
{
    double lambda = lambdaV(x, y);
    
    double a1 = lambda*h_y/h_x;
    double a2 = lambda*h_x/h_y;

    //верхний треугольник
    G[0][1] = -2 * a1 + a2;
    G[0][2] = a1 - 2 * a2;
    G[0][3] = -a1 - a2;
    G[1][2] = -a1 - a2;
    G[1][3] = a1 + 2 * a2;
    G[2][3] = -2 * a1 + a2;

    //диагональ
    G[0][0] = 2 * a1 + 2 * a2;
    G[1][1] = 2 * a1 + 2 * a2;
    G[2][2] = 2 * a1 + 2 * a2;
    G[3][3] = 2 * a1 + 2 * a2;

    //нижний треугольник
    G[1][0] = -2 * a1 + a2;
    G[2][0] = a1 - 2 * a2;
    G[2][1] = a1 - a2;
    G[3][0] = a1 - a2;
    G[3][1] = a1 - 2 * a2;
    G[3][2] = -2 * a1 + a2;
}

//доделать гамму
void LocalM_matrix() //построение локальной матрицы масс
{
    double gamma;
    double a = (gamma * h_x * h_y) / 36.;
    //верхний треугольник
    M[0][1] = 2 * a;
    M[0][2] = 2 * a;
    M[0][3] = 2 * a;
    M[1][2] = a;
    M[1][3] = 2 * a;
    M[2][3] = 2 * a;

    //диагональ
    M[0][0] = 4 * a;
    M[1][1] = 4 * a;
    M[2][2] = 4 * a;
    M[3][3] = 4 * a;

    //нижний треугольник
    M[1][0] = 2 * a;
    M[2][0] = 2 * a;
    M[2][1] = a;
    M[3][0] = a;
    M[3][1] = 2 * a;
    M[3][2] = 2 * a;
}

//построение локального вектора правых частей
void LocalB_vector(double f1, double f2, double f3, double f4) 
{
    double a = (h_x * h_y) / 36.;
    b[0] = a * (4 * f1 + 2 * f2 + 2 * f3 + f4);
    b[1] = a * (2 * f1 + 4 * f2 + f3 + 2 * f4);
    b[2] = a * (2 * f1 + f2 + 4 * f3 + 2 * f4);
    b[3] = a * (f1 + 2 * f2 + 2 * f3 + 4 * f4);
}
int main() {
    int N = 9; // Количество узлов
    int Kel = 1; // Количество элементов

    Portret_builders();

    // Вывод для проверки
    cout << "ig: ";
    for (int i = 0; i <= N; ++i) {
        cout << ig[i] << " ";
    }
    cout << endl;

    cout << "jg: ";
    for (size_t i = 0; i < jg.size(); ++i) {
        cout << jg[i] << " ";
    }
    cout << endl;

    return 0;
}