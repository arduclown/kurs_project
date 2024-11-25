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

int x_size = 0; //кол-во элементов по оси х
int y_size = 0; //кол-во элементов по оси у
int global_size = 0; //общее кол-во элементов
double h_x = 1; //шаг сетки по х
double h_y = 1; //шаг сетки по у
int L = 0; //кол-во подобластей
int cntKraev = 0; //кол-во краевых условий

vector<int> X_grid; //разбиение на сетку по х (координаты)
vector<int> Y_grid; //разбиение на сетку по у (координаты)
vector<vector<double> > G_l; //матрица жесткости
vector<vector<double> > M_l; //матрица масс
vector<double> b_l; //вектор правых частей
vector<vector<double> > localA; //матрица коэффициентов
vector<double> localB; //вектор коэффициентов
vector<double> B; //вектор глобальный
vector<Rectangle> rectangles;


//хранение матрицы в разреженном строчном формате
vector<double> di; //диагональ
vector<double> gg; //ненулевые компоненты матрицы
vector<int> ig; //индексы строк ненулевых компонет
vector<int> jg; //индексы столбцов ненулевых компонент

//краевые условия
vector<vector<int> > kraevye_uslov; //кравеые условия
vector<vector<double> > val_A3; //значения 3 краевых условий для матрицы
vector<double> val_B; //значения 3 и 2 краевых условий для вектора b

//решение слау
vector<double> x0; //начальное приближение
vector<double> y;
vector<double> z;
vector<double> t;
vector<double> r;
int maxiter = 10000;
double e = 1e-14;


double Func(double x, double y) 
{
    return 0;
}

double lambdaV ()
{
    return 1;
}

double gammaV()
{
    return 0;
}

double Betta(int num)
{
    return 1;
}

//доделать учет 
double u_Betta(int num)
{
    return 1;
}

double Tetta(int num)
{
    return 1;
}

double UG(int num)
{
    return 1;
}

/*получение глобального номера узла (из учебника)*/
int IndexOfUnknown(int i, int j) { //возвращает глобальный номер j-й локальной базисной функции внутри i-го конечного элемента
    // i (ielem) - номер конечного элемента
    // j (localIndex) - локальный индекс узла внутри элемента (0, 1, 2 или 3)

    //элементы нумеруются последовательно по строкам (слева направо, сверху вниз)
    int elementsPerRow = x_size-1; // Количество ребер в строке. ( - 1, т.к. каждое ребро соединяет два узла)
    int row = i / elementsPerRow; // Определение строки
    int col = i % elementsPerRow; // Определение столбца

    // Возвращаем глобальный индекс узла на основе его локального индекса
    switch (j) {
        case 0: //Левый нижний узел
            return row * x_size + col;
        case 1: //Правый нижний узел
            return row * x_size + col + 1;
        case 2: //Левый верхний узел
            return (row + 1) * x_size + col;
        case 3: //Правый верхний узел
            return (row + 1) * x_size + col + 1;
        default:
            cout << "Invalid local index for element" << endl;
            return -1; // если индекс выходит за пределы допустимого диапазона
    }

}

/*построение сетки*/
void GridBuilder()
{
    ifstream input("matrix.txt");
    input >> x_size >> y_size; 
    X_grid.resize(x_size);
    Y_grid.resize(y_size);
    for (int i = 0; i < x_size; i++)
        input >> X_grid[i];
    for (int i = 0; i < y_size; i++)
        input >> Y_grid[i];
    global_size = x_size * y_size;

    //построение подобластей    
    input >> L;
    rectangles.resize(L); //здесь лежит информация о подобластях.
    for (int i = 0; i < L; i++)
        input >> rectangles[i].position >> rectangles[i].x_left >> rectangles[i].x_right >> rectangles[i].y_down >> rectangles[i].y_up;
    input.close();
}

/*построение локальной матрицы жесткости*/
vector<vector<double> > LocalG_matrix() 
{
    vector<vector<double> > G(4, vector<double>(4));
    double lambda = 1;
    
    double a1 = lambda*h_y/(6.*h_x);
    double a2 = lambda*h_x/(6.*h_y);

    //верхний треугольник
    G[0][1] = -2 * a1 + a2;
    G[0][2] = a1 - 2 * a2;
    G[0][3] = -a1 - a2;
    G[1][2] = -a1 - a2;
    G[1][3] = a1 - 2 * a2;
    G[2][3] = -2 * a1 + a2;

    //диагональ
    G[0][0] = 2 * a1 + 2 * a2;
    G[1][1] = 2 * a1 + 2 * a2;
    G[2][2] = 2 * a1 + 2 * a2;
    G[3][3] = 2 * a1 + 2 * a2;

    //нижний треугольник
    G[1][0] = -2 * a1 + a2;
    G[2][0] = a1 - 2 * a2;
    G[2][1] = -a1 - a2;
    G[3][0] = -a1 - a2;
    G[3][1] = a1 - 2 * a2;
    G[3][2] = -2 * a1 + a2;

    return G;
}

/*построение локальной матрицы масс*/
vector<vector<double> > LocalM_matrix() 
{
    vector<vector<double> > M(4, vector<double>(4));
    double gamma = 2;
    double a = (gamma * h_x * h_y) / 36.;
    //верхний треугольник
    M[0][1] = 2 * a;
    M[0][2] = 2 * a;
    M[0][3] = a;
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
    
    return M;
}

/*построение локального вектора правых частей*/
vector<double> LocalB_vector(int num_L) 
{
    vector<double> b(4);
    double f1 = Func(rectangles[num_L].x_left, rectangles[num_L].y_down);
    double f2 = Func(rectangles[num_L].x_right, rectangles[num_L].y_down);
    double f3 = Func(rectangles[num_L].x_left, rectangles[num_L].y_up);
    double f4 = Func(rectangles[num_L].x_right, rectangles[num_L].y_up);
    double a = (h_x * h_y) / 36.;
    b[0] = a * (4 * f1 + 2 * f2 + 2 * f3 + f4);
    b[1] = a * (2 * f1 + 4 * f2 + f3 + 2 * f4);
    b[2] = a * (2 * f1 + f2 + 4 * f3 + 2 * f4);
    b[3] = a * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

    return b;
}

/*вычисление локальной матрицы*/
vector<vector<double> > LocalMatrix(int num_L) //номер конечного элемента
{
    vector<vector<double> > A (4, vector<double>(4));

    G_l = LocalG_matrix();
    M_l = LocalM_matrix();
    for (int i = 0; i < G_l.size(); i++)
    {
        for(int j = 0; j < M_l.size(); j++)
            A[i][j] = G_l[i][j] + M_l[i][j];
    }
    return A;

}

/*построение портрета глоабльной матрицы*/
void Portret_builders()
{

    vector<set<int> > list; // временный массив для хранения списка
    list.resize(global_size);

    for (int k = 0; k < L; k++)
    {
        vector<int> tmp; //тут будут хранится глоабльные индексы узлов
        for (int p = 0; p < 4; p++)
            tmp.push_back(IndexOfUnknown(k, p));
        for (int i = 0; i < 4; i++)
        {
            int ind1 = tmp[i];
            for (int j = i + 1; j < 4; j++)
            {
                int ind2 = tmp[j];
                list[ind1].insert(ind2);
                list[ind2].insert(ind1);
            }
        }
    }
    // Заполнение массива ig (индексов строк)
    ig.resize(global_size + 1);
    ig[0] = 0;
    ig[1] = 0;
    for (int i = 0; i < global_size; i++) {
        ig[i+1] = ig[i] + list[i].size();
    }

    // Заполнение массива jg (индексов столбцов)
    jg.resize(ig[global_size]);
    for (int i = 0, k = 0; i < global_size; i++) {
        for (int j : list[i]) {
            jg[k] = j;
            k++;
        }
    }

    di.resize(global_size);
    gg.resize(ig[global_size]);
}

/*построение глобальной матрицы*/
void GlobalMatrix(vector<vector<double> > &localA, int k) //индекс конечного элемента
{
    vector<int> tmp; //тут будут хранится глоабльные индексы узлов
        for (int p = 0; p < 4; p++)
        {
            tmp.push_back(IndexOfUnknown(k, p));
            di[tmp[p]] += localA[p][p];
        }
        // Добавление вненедиагональных элементов
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < i; j++) {
            // Глобальные индексы узлов
            int igg = tmp[i];
            int jgg = tmp[j];

            // Обеспечение порядка (igg > jgg)
            if (igg < jgg) {
                swap(igg, jgg);
            }

            // Поиск позиции в массиве разреженной матрицы
            int index = ig[igg];

            while(index < ig[igg + 1] && jg[index]!=jgg) index++; 
    

            // Добавление элементов в разреженные массивы
            gg[index] += localA[i][j]; 
        }
    }
}

/*сборка глобального вектора*/
vector<double> BuildGlobalVector(vector<double> &l_B, int k) {

    l_B=LocalB_vector(k);
    for (int i = 0; i < 4; i++) {
        int ind = IndexOfUnknown(k, i);
        B[ind] += l_B[i];
    }
    return B;
}

/*построение краевых условий*/
void Kraevye() //описание: 1 - тип кравевого усл, 2 - номер границы, 3 - координата x начало, 4 - координата x конец, 5 - координата y начало, 6 - координата y конец, 
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
            //cout << kraevye_uslov[i][j] << " ";
        }
        //cout << endl;
    }
    input.close();
}


void Third_K(int k, const vector<int> &kraevye_uslov) {
    vector<int> tmp; // Узлы на границе

    // Коэффициенты для третьего краевого условия
    double betta = Betta(kraevye_uslov[1]);
    double u_betta1 = u_Betta(kraevye_uslov[1]);
    double u_betta2 = u_Betta(kraevye_uslov[1]);

    // Координаты границы
    double x_left = kraevye_uslov[2];
    double x_right = kraevye_uslov[3];
    double y_down = kraevye_uslov[4];
    double y_up = kraevye_uslov[5];

    // Проверка узлов текущего конечного элемента
    for (int localNode = 0; localNode < 4; localNode++) {
        int globalNode = IndexOfUnknown(k, localNode);

        // Получение координат узла из X_grid и Y_grid
        double x_coord = X_grid[globalNode % x_size]; // X-координата
        double y_coord = Y_grid[globalNode / x_size]; // Y-координата

        // Проверяем, попадает ли узел на границу
        if ((x_coord == x_left && x_coord == x_right) || 
            (y_coord == y_down && y_coord == y_up)) {
            tmp.push_back(globalNode); // Добавляем узел на границе
        }
    }

    // Проверка: граница должна быть между двумя узлами
    if (tmp.size() != 2) {
        cerr << "Ошибка: граница не соответствует двум узлам!" << endl;
        return;
    }

    // Длина отрезка
    double h;
    double x1 = X_grid[tmp[0] % x_size];
    double x2 = X_grid[tmp[1] % x_size];
    double y1 = Y_grid[tmp[0] / x_size];
    double y2 = Y_grid[tmp[1] / x_size];
    h = (x1 == x2) ? abs(y2 - y1) : abs(x2 - x1); // Если x совпадает, длина по y, иначе по x

    // Коэффициент для третьего рода
    double a = (betta * h) / 6.0;

    val_A3.resize(4, vector<double>(4));
    val_A3[0][0] = 2 * a;
    val_A3[0][1] = a;
    val_A3[1][0] = a;
    val_A3[1][1] = 2 * a;

    val_B.resize(2);
    val_B[0] = a * (2 * u_betta1 + u_betta2);
    val_B[1] = a * (u_betta1 + 2 * u_betta2);

    // Занесение в глобальную матрицу и вектор
    for (int i = 0; i < 2; i++) {
        int global_i = tmp[i];

        // Обновление диагональных элементов
        di[global_i] += val_A3[i][i];

        // Обновление внедиагональных элементов
        for (int j = 0; j < 2; j++) {
            if (i != j) {
                int global_j = tmp[j];

                // Найти индекс в разреженной структуре
                int index = ig[global_i];
                while (index < ig[global_i + 1] && jg[index] != global_j) {
                    index++;
                }

                // Добавить значение в разреженную структуру
                if (jg[index] == global_j) {
                    gg[index] += val_A3[i][j];
                }
            }
        }

        // Обновление правой части
        b_l[global_i] += val_B[i];
    }

    
}


/*учёт 2х краевых условий*/
void Second_K(int k, const vector<int> &kraevye_uslov)
{
    vector<int> tmp; // Узлы на границе
    double h; // Длина отрезка (либо по x, либо по y)

    // Значения потока на границе
    double tetta1 = Tetta(kraevye_uslov[1]); // Поток в первом узле
    double tetta2 = Tetta(kraevye_uslov[1]); // Поток во втором узле

    // Координаты границы
    double x_left = kraevye_uslov[2];
    double x_right = kraevye_uslov[3];
    double y_down = kraevye_uslov[4];
    double y_up = kraevye_uslov[5];

    // Проверка узлов текущего конечного элемента
    for (int localNode = 0; localNode < 4; localNode++) {
        int globalNode = IndexOfUnknown(k, localNode);

        // Получение координат узла из X_grid и Y_grid
        double x_coord = X_grid[globalNode % x_size]; // X-координата
        double y_coord = Y_grid[globalNode / x_size]; // Y-координата

        // Проверяем, попадает ли узел на границу
        if ((x_coord == x_left && x_coord == x_right) || 
            (y_coord == y_down && y_coord == y_up)) {
            tmp.push_back(globalNode); // Добавляем узел на границе
        }
    }

    // Проверка: граница должна быть между двумя узлами
    if (tmp.size() != 2) {
        cerr << "Ошибка: граница не соответствует двум узлам!" << endl;
        return;
    }

    // Длина отрезка
    double x1 = X_grid[tmp[0] % x_size];
    double x2 = X_grid[tmp[1] % x_size];
    double y1 = Y_grid[tmp[0] / x_size];
    double y2 = Y_grid[tmp[1] / x_size];
    h = (x1 == x2) ? abs(y2 - y1) : abs(x2 - x1); // Если x совпадает, длина по y, иначе по x

    // Рассчёт локального вектора нагрузки
    double a = h / 6.; // Коэффициент для второго рода
    val_B.resize(2);
    val_B[0] = a * (2 * tetta1 + tetta2);
    val_B[1] = a * (tetta1 + 2 * tetta2);

    b_l[tmp[0]] += val_B[0];
    b_l[tmp[1]] += val_B[1];
}

/*учёт 1х краевых условий*/
void First_K(int k, const vector<int> &kraevye_uslov) {
    vector<int> tmp; // Узлы на границе
    double u_g = UG(kraevye_uslov[1]); // Значение на границе

    
    // Координаты границы
    double x_left = kraevye_uslov[2];
    double x_right = kraevye_uslov[3];
    double y_down = kraevye_uslov[4];
    double y_up = kraevye_uslov[5];

    // Проверка узлов текущего конечного элемента
    for (int localNode = 0; localNode < 4; localNode++) {
        int globalNode = IndexOfUnknown(k, localNode);

        // Получение координат узла из X_grid и Y_grid
        double x_coord = X_grid[globalNode % x_size]; // X-координата
        double y_coord = Y_grid[globalNode / x_size]; // Y-координата

        // Проверяем, попадает ли узел на границу
        if ((x_coord == x_left && x_coord == x_right) || 
            (y_coord == y_down && y_coord == y_up) ) {
            tmp.push_back(globalNode); // Добавляем узел на границе
        }
        cout << "Глобальный узел: " << globalNode << " Координаты: (" << x_coord << ", " << y_coord << ")" << endl;
        cout << "Граница: x = [" << x_left << ", " << x_right << "], y = [" << y_down << ", " << y_up << "]" << endl;

    }
    
    // Обработка узлов на границе
    for (int globalNode : tmp) {
        // Установка диагонального элемента
        di[globalNode] = 1;

        // Установка правой части
        B[globalNode] = u_g;

        // Зануление строки
        for (int i = ig[globalNode]; i < ig[globalNode + 1]; i++) {
            gg[i] = 0; // Обнуление всех недиагональных элементов строки
        }

        // Зануление столбца
        for (int i = 0; i < global_size; i++) {
            for (int j = ig[i]; j < ig[i + 1]; j++) {
                if (jg[j] == globalNode) {
                    gg[j] = 0; // Обнуление столбца
                }
            }
        }
    }
}


void AddKraevye(int num, int k, vector<int> &kraevye_uslov) 
{
    //num - тип краевого условия
    //k - номер конечного элемента
    switch(num){
        case(1):
            First_K(k, kraevye_uslov);
            break;
        case(2):
            Second_K(k, kraevye_uslov);
            break;
        case(3):
            Third_K(k, kraevye_uslov);
            break;
        default:
            cout << "Incorrect Boundary Condition Number" << endl;
            break;
    }
}

/*сборка глоабльной матрицы!!!:3*/
void BuildGlobalMatrix() {
    Portret_builders();
    // Проход по всем конечным элементам
    for (int k = 0; k < L; k++) {
        localA=LocalMatrix(k);
        GlobalMatrix(localA, k);
        //b_l=LocalB_vector(k);
        B = BuildGlobalVector(b_l, k);
        for (int i = 0; i < cntKraev; i++)
            AddKraevye(kraevye_uslov[i][0], k, kraevye_uslov[i]);
    }
}

/*вывод глобальной матрицы*/
void PrintGlobalMatrix(int global_size) {
    // Создаем временную двумерную матрицу для удобного отображения
    vector<vector<double> > fullMatrix;
    fullMatrix.resize(global_size, vector<double>(global_size, 0.0));

    // Заполнение диагональных элементов
    for (int i = 0; i < global_size; i++) {
        fullMatrix[i][i] = di[i];
    }

    // Заполнение недиагональных элементов
    for (int i = 0; i < global_size; i++) {
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            int col = jg[j];        // Столбец из разреженной структуры
            fullMatrix[i][col] = gg[j]; // Значение из массива gg
            fullMatrix[col][i] = gg[j]; // Симметричность матрицы
        }
    }

    // Вывод матрицы в консоль
    cout << "Global Matrix:" << endl;
    for (int i = 0; i < global_size; i++) {
        for (int j = 0; j < global_size; j++) {
            cout << setw(10) << fullMatrix[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Global Vector:" << endl;
    for (int i = 0; i < global_size; i++) {
        cout << setw(10) << B[i] << " ";
    }
    cout << endl;
}
double Norma(vector<double> &v)
{
	double norma = 0;
	for (int i = 0; i < global_size; i++)	
		norma += v[i] * v[i];
	return sqrt(norma);
}

double vector_multiplication(vector<double> &v1, vector<double> &v2)
{
	double sum = 0;
	for (int i = 0; i < global_size; i++)
		sum += v1[i] * v2[i];
	return sum;
}

vector<double> matrix_on_vector_multiplication(vector<double> &v1, vector<double> &v2) 
{
	for (int i = 0; i < global_size; i++)
		v2[i] = di[i] * v1[i];
	for (int i = 1; i < global_size; i++)
	{
		int i0 = ig[i];
		int i1 = ig[i + 1];
		for (int j = 0; j < (i1 - i0); j++)
		{
			v2[i] += gg[i0+j] * v1[jg[i0+j]];
			v2[jg[i0+j]] += gg[i0+j] * v1[i];
		}
		
	}
	return v2;
}

vector<double> vector_sum(vector<double> &v1, vector<double> &v2, vector<double> &v3, double a)
{
	for (int i = 0; i < global_size; i++)
		v3[i] = v1[i] + a * v2[i];
	return v3;
}

void revDiagonal(vector<double>& v1, vector<double>& v2)	//	b = a/di;
{
	for (int i = 0; i < global_size; i++)
		v1[i] = v2[i] / di[i];
}

// МСГ с предобусловливанием
void MSG_P()
{
	double alpha, betta;
	double scal_rk;
	int k_iter = 0;
    for (int i = 0; i < global_size; i++) x0[i] = 0;
	y = matrix_on_vector_multiplication(x0, y); // y = Ax^0
	vector_sum(B, y, r, -1);	//	r0 = f - Ax0
	//пусть M = D, где D - диагональ матрицы A
	revDiagonal(z, r); // z0 = r0 * (1 / D), где 1/D = M^-1
	for (int k = 1; k < maxiter; k++)
	{
		revDiagonal(t, r); // M^-1*r^(k-1)
		y = matrix_on_vector_multiplication(z, y); // y = A*z_(k-1)
		scal_rk = vector_multiplication(t, r); //(M^-1*r^(k-1),r^(k-1))
		alpha = scal_rk / vector_multiplication(y, z); // alpha_k = (M^-1*r^(k-1),r^(k-1))) / (A*z_(k-1),z_(k-1))
		x0 =vector_sum(x0, z, x0, alpha); // x_k = x_(k-1) + alpha_k * z_(k-1)
		y = matrix_on_vector_multiplication(z, y); // A*z_(k-1) = y 
		vector_sum(r, y, r, -alpha); // r_k = r_(k-1) - alpha_k * A * z_(k-1)
		revDiagonal(t, r); // t = M^-1*r^k
		betta = vector_multiplication(t, r) / scal_rk; //	betta_k = (r_k / D,r_k) / (r_(k-1),r_(k-1))
		vector_sum(t, z, z, betta);//	z_k = M^-1*r_k + betta_k * z_(k-1)
		//betta = 0.;
		if (Norma(r) / Norma(B) < e) // ||r_k|| / ||f|| < e
		{
			cout << Norma(r) << " " << Norma(B) << endl;
			k_iter = k;
			cout << k_iter;
			return;
		}
	}
	k_iter = maxiter;
	cout << k_iter;
}

void Output()
{
	//cout << fixed << setprecision(15);
    //cout.precision(15); 
	for (int i = 0; i < global_size; i++)
		cout << x0[i] << endl;
}

void SLAU()
{
    x0.resize(global_size);
    y.resize(global_size);
    z.resize(global_size);
    t.resize(global_size);
    r.resize(global_size);
    MSG_P();
    Output();
}

int main() {
    GridBuilder();
    G_l.resize(16);
    M_l.resize(16);
    localA.resize(16);
    b_l.resize(4);
    B.resize(global_size);
    Kraevye();
    //cout << cntKraev << endl;
    BuildGlobalMatrix();
    //cout << "we here" << endl;
    SLAU();
    PrintGlobalMatrix(global_size);

    return 0;
}