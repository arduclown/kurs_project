#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
using namespace std;

struct Rectangle
{
    int position; //подобласть 
    int x_left; 
    int x_right;  //координаты кравеого
    int y_down;
    int y_up;
};

int x_size = 0; //кол-во элементов по оси х
int y_size = 0; //кол-во элементов по оси у
int global_size = 0; //общее кол-во элементов
double h_x = 0; //шаг сетки по х
double h_y = 0; //шаг сетки по у
int L = 0; //кол-во подобластей
int cntKraev = 0; //кол-во краевых условий

vector<int> X_grid; //разбиение на сетку по х (координаты)
vector<int> Y_grid; //разбиение на сетку по у (координаты)
vector<vector <int> > GlobalNodes;
vector<vector<double> > G_l; //матрица жесткости
vector<vector<double> > M_l; //матрица масс
vector<double> b_l; //вектор правых частей
vector<vector<double> > localA; //матрица коэффициентов
vector<double> localB; //вектор коэффициентов
vector<double> B; //вектор глобальный
vector<Rectangle> rectangles;

//для разбиений
vector<vector<double> > nodes;
vector<double> X;
vector<double> Y;
vector<int> IX;
vector<int> IY;
int Nx, Ny;

//хранение матрицы в разреженном строчном формате
vector<double> di; //диагональ
//vector<double> gg; //ненулевые компоненты матрицы
vector<double> al;
vector<double> au;
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
vector<double> p;
vector<double> x;

int maxiter = 10000;
double e = 1e-16;


double Func(int num, double x, double y) 
{
   switch(num){
        case 1: return 10 * x * x * x - 60;
        case 2: return 1.8 + 0.1 * x;
        case 3: 0.0;
   }
}

double lambdaV (int num)
{
    switch(num){
        case 1: return 10;
        case 2: return 10;
        case 3: return 10;
    }
}

double gammaV(int num)
{
    switch(num){
        case 1: return 10;
        case 2: return 10;
        case 3: return 10;
    }

}

double BettaV(int num)
{
    switch(num)
    {
        case 1:
            return 1;
            //return x;
        case 2:
            return 1;
            //return y + 4;
        case 3:
            return 0.5;
    }
}

double u_Betta(int num, double x, double y)
{
    //num - номер границы
    switch(num)
    {
        case 1:
            return x * x * x;
            //return x;
        case 2:
            return y-9;
            //return y + 4;
        case 3:
            return -1;
    }
}

double Tetta(int num, double x, double y)
{
    //num - номер границы
    switch(num)
    {
        case 1:
            return 30 * x * x;
        case 2:
            return -10;
    }
}

double UG(int num, double x, double y)
{
    //num - номер границы
    switch(num)
    {
        case 1:
            return x * x * x;
        case 2:
            return x * x * x;
        case 3:
            return 0.1*x + 1.8;
    }


}


/*получение глобального номера узла (из учебника)*/
int IndexOfUnknown(int i, int j) {
    // i (ielem) - номер конечного элемента
    // j (localIndex) - локальный индекс узла внутри элемента (0, 1, 2 или 3)
    
    int elementsPerRow = x_size - 1; // Количество элементов в строке (Nx - 1)

    int row = i / elementsPerRow;   // Определение строки
    int col = i % elementsPerRow;   // Определение столбца

    // Глобальный индекс узла на основе локального индекса
    int globalIndex;
    switch (j) {
        case 0: // Левый нижний узел
            globalIndex = row * x_size + col;
            break;
        case 1: // Правый нижний узел
            globalIndex = row * x_size + col + 1;
            break;
        case 2: // Левый верхний узел
            globalIndex = (row + 1) * x_size+ col;
            break;
        case 3: // Правый верхний узел
            globalIndex = (row + 1) * x_size+ col + 1;
            break;
        default:
            cout << "Invalid local index for element" << endl;
            return -1; // Ошибка, если индекс выходит за пределы допустимого диапазона
    }

    return globalIndex;
}

/*построение сетки*/
void GridBuilder()
{
    ifstream input("./test/test5/matrix.txt");
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

    nodes.resize(global_size);
    int i = 0;

    for (double y_c : Y_grid) {
        for (double x_c : X_grid) {
            nodes[i].push_back(x_c);
            nodes[i].push_back(y_c); // Добавляем узел (x, y)
            i++;
        }
    }
}

/*построение локальной матрицы жесткости*/
vector<vector<double> > LocalG_matrix(int k, int num_L) 
{
    vector<vector<double> > G(4, vector<double>(4));
    double lambda = lambdaV(k);

    h_x = rectangles[num_L].x_right- rectangles[num_L].x_left;
    h_y = rectangles[num_L].y_up - rectangles[num_L].y_down;
    double a1 = (lambda/6)*(h_y/h_x);
    double a2 = lambda*h_x/(6*h_y);

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
vector<vector<double> > LocalM_matrix(int k, int num_L) 
{
    vector<vector<double> > M(4, vector<double>(4));
    h_x = rectangles[num_L].x_right- rectangles[num_L].x_left;
    h_y = rectangles[num_L].y_up - rectangles[num_L].y_down;
    double gamma = gammaV(k);
    double a = (gamma * h_x * h_y) / 36;
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
    //num_L - номер краевого
    h_x = rectangles[num_L].x_right- rectangles[num_L].x_left;
    h_y = rectangles[num_L].y_up - rectangles[num_L].y_down;
    vector<double> b(4);
    double f1 = Func(rectangles[num_L].position, rectangles[num_L].x_left, rectangles[num_L].y_down);
    double f2 = Func(rectangles[num_L].position, rectangles[num_L].x_right, rectangles[num_L].y_down);
    double f3 = Func(rectangles[num_L].position,rectangles[num_L].x_left, rectangles[num_L].y_up);
    double f4 = Func(rectangles[num_L].position,rectangles[num_L].x_right, rectangles[num_L].y_up);
    double a = (h_x * h_y) / 36;
    b[0] = a * (4 * f1 + 2 * f2 + 2 * f3 + f4);
    b[1] = a * (2 * f1 + 4 * f2 + f3 + 2 * f4);
    b[2] = a * (2 * f1 + f2 + 4 * f3 + 2 * f4);
    b[3] = a * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

    
    cout << "VECTOR number of elem: " << num_L + 1 << endl;
    for (int i = 0; i < 4; i++) cout << b[i] << " ";
    cout << endl;
    cout << endl;
    return b;
}

/*вычисление локальной матрицы*/
vector<vector<double> > LocalMatrix(int num_L, int k) 
{
    //numL - номер подобласти
    vector<vector<double> > A (4, vector<double>(4));

    G_l = LocalG_matrix(num_L, k);
    M_l = LocalM_matrix(num_L, k);
    cout << "number elem: " << k+1 << endl;
    for (int i = 0; i < G_l.size(); i++)
    {
        for(int j = 0; j < M_l.size(); j++){
            A[i][j] = G_l[i][j] + M_l[i][j];
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    return A;

}

/*построение портрета глоабльной матрицы*/
void Portret_builders()
{

    vector<set<int> > list; // временный массив для хранения списка связности
    list.resize(global_size);

    int m = 0;
    for (; m < L; m++)
    {
        cout << m << endl;
        vector<int> tmp; //тут будут хранится глоабльные индексы узлов
        for (int p = 3; p >= 0; p--)
            //int iN = IndexOfUnknown(m,p);
            //нужно реализовать проверку на фиктивный узел. Если фиктивный, то не выдавать 
            //if (!IsFictionNode(m, ))
            tmp.push_back(IndexOfUnknown(m, p));
        for (int i = 0; i < 4; i++)
        {
            int ind1 = tmp[i];
            for (int j = i + 1; j < 4; j++)
            {
                int ind2 = tmp[j];
                list[ind1].insert(ind2);
                //list[ind2].insert(ind1);
            }
        }
    }
    
    // Заполнение массива ig (индексов строк)
    ig.resize(global_size + 1);
    ig[0] = 0;
    ig[1] = 0;
    for (int i = 0; i < global_size; i++)
    {
        ig[i + 1] = ig[i] + list[i].size();
    }

    // Заполнение массива jg (индексов столбцов)
    jg.resize(ig[global_size]);
    int k = 0;              // Индекс для заполнения jg
    jg[k] = 0;
    for (int i = 0; i < global_size; i++) {
        for (int elem : list[i]) { // Используем итерацию по элементам множества
            jg[k] = elem;       // Добавляем элемент в jg
            k++;
        }
    }

    di.resize(global_size);
    al.resize(ig[global_size]);
    au.resize(ig[global_size]);
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
            al[index] += localA[i][j]; 
            au[index] += localA[i][j];
        }
    }
}

/*сборка глобального вектора*/
vector<double> BuildGlobalVector(vector<double> &l_B, int k) 
{
    //k - номер конечного элемента
    vector<double> tmp;
    tmp.resize(4);
    for(int i = 0; i < 4; i++) tmp[i] = IndexOfUnknown(k,i);
    for (int i = 0; i < 4; i++) {
        int ind = tmp[i];
        B[ind] += l_B[i];
        //cout << B[ind] << endl;
    }
    return B;
}

/*построение краевых условий*/
void Kraevye() //описание: 1 - тип кравевого усл, 2 - номер границы, 3 - узел 1, 4 - узел 2
{
    ifstream input("./test/test5/kraevie.txt");
    input >> cntKraev;
    kraevye_uslov.resize(cntKraev);
    for (int i = 0; i < cntKraev; i++)
    {
        kraevye_uslov[i].resize(4);
        for (int j = 0; j < 4; j++)
            input >> kraevye_uslov[i][j];
    }
    input.close();
}

/*учёт 3х краевых условий*/
void Third_K(const vector<int> &kraevye_uslov) {
    vector<int> tmp; // Узлы на границе

    int i = kraevye_uslov[2];
    int j = kraevye_uslov[3];

    // Коэффициенты для третьего краевого условия
    double betta = BettaV(kraevye_uslov[1]);
    double u_betta1 = u_Betta(kraevye_uslov[1], nodes[i][0], nodes[i][1]);
    double u_betta2 = u_Betta(kraevye_uslov[1], nodes[j][0], nodes[j][1]);

    tmp.push_back(i);
    tmp.push_back(j);

    // Длина отрезка
    double h;
    double x1 = nodes[i][0];
    double x2 = nodes[j][0];
    double y1 = nodes[i][1];
    double y2 = nodes[j][1];
    h = (x1 == x2) ? abs(y2 - y1) : abs(x2 - x1); // Если x совпадает, длина по y, иначе по x

    // Коэффициент для третьего рода
    double a = (betta * h) / 6.0;

    val_A3.resize(2, vector<double>(2));
    val_A3[0][0] = 2 * a;
    val_A3[0][1] = a;
    val_A3[1][0] = a;
    val_A3[1][1] = 2 * a;

    val_B.resize(2);
    val_B[0] = a * (2 * u_betta1 + u_betta2);
    val_B[1] = a * (u_betta1 + 2 * u_betta2);

    B[i] += val_B[0];
    B[j] += val_B[1];

    di[i] += val_A3[0][0];
    di[j] += val_A3[1][1];

    int ia, ja;
    ia = j;
    ja = i;
    int index = ig[ia];
    int flag = 1;
    for (; index < ig[ia + 1] && flag; index++)
		if (jg[index] == ja) flag = 0;
	index--;
	al[index] += val_A3[1][0];
	au[index] += val_A3[0][1];
}

/*учёт 2х краевых условий*/
void Second_K(const vector<int> &kraevye_uslov)
{
    vector<int> tmp; // Узлы на границе
    double h; // Длина отрезка (либо по x, либо по y)

    
    int i = kraevye_uslov[2];
    int j = kraevye_uslov[3];

    // Значения потока на границе
    double tetta1 = Tetta(kraevye_uslov[1], nodes[i][0], nodes[i][1]); // Поток в первом узле
    double tetta2 = Tetta(kraevye_uslov[1], nodes[j][0], nodes[j][1]);

    tmp.push_back(i);
    tmp.push_back(j);

    // Длина отрезка
    double x1 = nodes[i][0];
    double x2 = nodes[j][0];
    double y1 = nodes[i][1];
    double y2 = nodes[j][1];
    h = (x1 == x2) ? abs(y2 - y1) : abs(x2 - x1); // Если x совпадает, длина по y, иначе по x

    // Рассчёт локального вектора нагрузки
    double a = h / 6.0; // Коэффициент для второго рода
    val_B.resize(2);
    val_B[0] = a * (2 * tetta1 + tetta2);
    val_B[1] = a * (tetta1 + 2 * tetta2);

    B[i] += val_B[0];
    B[j] += val_B[1];
}

/*учёт 1х краевых условий*/
void First_K(const vector<int> &kraevye_uslov) {
    vector<int> tmp; // Узлы на границе
    int p = kraevye_uslov[2];
    int s = kraevye_uslov[3];
    tmp.push_back(p);
    tmp.push_back(s);

    // Обработка узлов на границе
    for (int globalNode : tmp) {
        double x_coord = nodes[globalNode][0]; // X-координата
        double y_coord = nodes[globalNode][1]; // Y-координата
        double u_g = UG(kraevye_uslov[1], x_coord, y_coord);
        // диагональ
        di[globalNode] = 1;

        // правая часть
        B[globalNode] = u_g;

        // зануление строки
        for (int i = ig[globalNode]; i < ig[globalNode + 1]; i++) {
            al[i] = 0; // Обнуление всех недиагональных элементов строки
        }
        
        for (int i = ig[globalNode]; i < ig[global_size]; i++)
        {
            if (jg[i] == globalNode) 
                au[i] = 0;
                
        }
    }
}

/*вывод глобальной матрицы*/
void PrintGlobalMatrix(int global_size) {
    // Создаем временную двумерную матрицу для удобного отображения
    vector<vector<double> > fullMatrix;
    fullMatrix.resize(global_size, vector<double>(global_size, 0));

    // Заполнение диагональных элементов
    for (int i = 0; i < global_size-1; i++) {
        fullMatrix[i][i] = di[i];
    }

    // Заполнение недиагональных элементов
    for (int i = 0; i < global_size-1; i++) {
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            int col = jg[j];        // Столбец из разреженной структуры
            fullMatrix[i][col] = al[j]; // Значение из массива gg
            fullMatrix[col][i] = au[j]; // Симметричность матрицы
        }
    }

    // Вывод матрицы в консоль
    cout << "Global Matrix:" << endl;
    for (int i = 0; i < global_size-1; i++) {
        for (int j = 0; j < global_size-1; j++) {
            cout << setw(10) << fullMatrix[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Global Vector:" << endl;
    for (int i = 0; i < global_size-1; i++) {
        cout << setw(10) << B[i] << " ";
    }
    cout << endl;
}

void AddKraevye(int num, vector<int> &kraevye_uslov) 
{
    //num - тип краевого условия
    //k - номер конечного элемента
    switch(num){
        case(1):
            First_K(kraevye_uslov);
            break;
        case(2):
            Second_K(kraevye_uslov);
            break;
        case(3):
            Third_K(kraevye_uslov);
            break;
        default:
            cout << "Incorrect Boundary Condition Number" << endl;
            break;
    }
}

/*сборка глоабльной матрицы!!!:3*/
void BuildGlobalMatrix() {
    Portret_builders();
    G_l.resize(16);
    M_l.resize(16);
    localA.resize(16);
    b_l.resize(4);
    // Проход по всем конечным элементам
    for (int k = 0; k < L; k++) {
        localA=LocalMatrix(rectangles[k].position, k);
        GlobalMatrix(localA, k);
        b_l=LocalB_vector(k);
        B = BuildGlobalVector(b_l, k);
        PrintGlobalMatrix(global_size);
    }
    for (int i = 0; i < cntKraev; i++)
            AddKraevye(kraevye_uslov[i][0], kraevye_uslov[i]);
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
	for (int i = 0; i < global_size; i++)
	{
		for (int j = ig[i]; j < ig[i+1]; j++)
		{
			v2[i] += al[j] * v1[jg[j]];
			v2[jg[j]] += au[j] * v1[i];
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

double Norma(vector<double> &v)
{
	double norma = 0;
	for (int i = 0; i < global_size; i++)	
		norma += v[i] * v[i];
	return sqrt(norma);
}

void Calc_Zk(double a, vector<double> &v1)
{
	for (int i = 0; i < global_size; i++)
		z[i] = v1[i] + a * z[i];
}

void Calc_Pk(double a, vector<double> &v1)
{
	for (int i = 0; i < global_size; i++)
		p[i] = v1[i] + a * p[i];
}

void Calc_Xk_Rk_L(double alphaK)
{
	for (int i = 0; i < global_size; i++)
	{
		x[i] = x[i] + alphaK * z[i];
		r[i] = r[i] - alphaK * p[i];
	}
}

double calcResidual() {
    double normb = 0;
	for (int i = 0; i < global_size; i++) {
        normb += B[i] * B[i];
        r[i] = B[i] - di[i] * x[i];
        int m = ig[i + 1];
        for (int k = ig[i]; k < m; k++) {
            int j = jg[k];
            r[i] -= al[k] * x[j];
            r[j] -= au[k] * x[i];
        }
    }

    double mul = vector_multiplication(r, r);
    return sqrt(mul/normb);
}

/*Локально оптимальная схема с диагональным предобуславливанием*/
void LOS_Dd()
{
    double alpha, betta;
	double normPr = Norma(B);
	/*r0 = f - Ax0*/
	matrix_on_vector_multiplication(x0, y); // y = Ax^0
	//r0 = f - Ax0	
	for (int i = 0; i < global_size; i++) r[i] = B[i] - y[i];

	/*z0 = r0*/
	for (int i = 0; i < global_size; i++)	
		z[i] = r[i];

	for (int i = 0; i < global_size; i++)
		x[i] = x0[i];
	
	/*p0 = A*z0*/
	matrix_on_vector_multiplication(z, p); 

	int k = 0;
	//type scal_r = vector_multiplication(r,r);
	double discrepancy = calcResidual();

	for (; k < maxiter && discrepancy > e; k++)
	{
		/*alpha_k = (r_k-1,p_k-1) / (p_k-1,p_k-1)*/
		double scal_pk_rk = vector_multiplication(r, p); // (r^(k-1),p^(k-1))
		double scal_pk = vector_multiplication(p, p); // (p^(k-1),p^(k-1))
		alpha = scal_pk_rk / scal_pk; // alpha_k = (r_k-1,p_k-1) / (p_k-1,p_k-1)

		Calc_Xk_Rk_L(alpha);
		/*betta_k = -(r_k,Ar_k) / (p_k,p_k)*/
		matrix_on_vector_multiplication(r, y); // y =A*r_(k)
		double Ar_p = vector_multiplication(y, p); // (A*r_k,p_k)
		betta = -Ar_p / scal_pk; //	betta_k = (r_k,Ar_k) / (p_k,p_k)

		/*z_k = r_k + betta_k * z_k*/
		Calc_Zk(betta, r);
		/*p_k = p_k + betta_k * p_k*/
		Calc_Pk(betta, y);

		discrepancy = sqrt(vector_multiplication(r,r)/vector_multiplication(B,B));
		
	}
	discrepancy = calcResidual();
	cout << "Iteration: "<< k << " RelDiscrepancy: " << discrepancy << endl;
}



void Output(ofstream& out)
{
	out << fixed << setprecision(15);
    out.precision(15); 
	for (int i = 0; i < global_size; i++)
		out << x[i] << endl;
}

void SLAU()
{
    x0.resize(global_size);
    x.resize(global_size);
    y.resize(global_size);
    z.resize(global_size);
    t.resize(global_size);
    r.resize(global_size);
    p.resize(global_size);
    LOS_Dd();
}


int main() {
    GridBuilder();
    B.resize(global_size);
    Kraevye();
    BuildGlobalMatrix();
    PrintGlobalMatrix(global_size);
    cout << "we here" << endl;
    SLAU();
    ofstream out("out.txt");
    Output(out);
    return 0;
}