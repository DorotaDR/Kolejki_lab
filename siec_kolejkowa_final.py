import numpy as np
from geneticalgorithm import geneticalgorithm as ga

# import warnings
#
#
# #disable warnings
# warnings.filterwarnings("ignore")


default_lambdas_we = [1.8, 1.2, 1]

default_mi_ir = np.array([[4, 4, 4],
                          [1, 0.8, 1],
                          [0.25, 0.2, 0.3],
                          [3, 3 / 2, 2],
                          [3 / 2, 1, 0],
                          [0, 0, 2],
                          [2, 2, 2]])


class AirportNetwork:

    def __init__(self, lambda_we_r=default_lambdas_we, mi_ir=default_mi_ir):

        self.lambda_we_r = lambda_we_r

        self.mi_ir = mi_ir

        self._columns = ["Wejście",
                         "Lądowanie",
                         "Przegląd",
                         "Serwis",
                         "Tankowanie",
                         "Pasażerowie",
                         "Wylot",
                         "Wyjscie"]

        self.p_shortdist = np.zeros((len(self._columns), len(self._columns)))
        self.p_shortdist[0][1] = 1
        self.p_shortdist[1][2] = 1
        self.p_shortdist[2][3] = 0.1
        self.p_shortdist[2][4] = 0.2
        self.p_shortdist[2][5] = 0.7
        self.p_shortdist[3][2] = 1
        self.p_shortdist[4][5] = 1
        self.p_shortdist[5][6] = 1
        self.p_shortdist[6][7] = 1

        self.p_longdist = np.zeros((len(self._columns), len(self._columns)))
        self.p_longdist[0][1] = 1
        self.p_longdist[1][2] = 1
        self.p_longdist[2][3] = 0.1
        self.p_longdist[2][4] = 0.9
        self.p_longdist[3][2] = 1
        self.p_longdist[4][5] = 1
        self.p_longdist[5][6] = 1
        self.p_longdist[6][7] = 1

        self.p_transport = np.zeros((len(self._columns), len(self._columns)))
        self.p_transport[0][1] = 1
        self.p_transport[1][2] = 1
        self.p_transport[2][3] = 0.1
        self.p_transport[2][4] = 0.7
        self.p_transport[2][5] = 0.2
        self.p_transport[3][2] = 1
        self.p_transport[4][5] = 1
        self.p_transport[5][6] = 1
        self.p_transport[6][7] = 1

        # C1_ij - macierz kosztów obsługi klasy w systemie
        self.C1_ij = [[0.9, 0.9, 0.9],
                      [0.2, 0.1, 0.2],
                      [0.1, 0.1, 0.1],
                      [0.3, 0.4, 0.3],
                      [0.7, 0.6, 0.0],
                      [0, 0, 0.7],
                      [0.5, 0.5, 0.5]]

        # C2_i - macierz kosztów liczby niezajętych kanałów
        self.C2_i = [0.1, 0.4, 0.1, 0.9, 0.6, 0.4, 0.5]

        self.ro_ir = None
        self.K_ir = None
        self.Q_ir = None
        self.m_nzi = None


    @staticmethod
    def _calculate_lambdas_ir(p_table, lambda_we):
        """
        :param p_table: macierz przejść pewnej klasy
        :param lambda_we: strumienie wejściowe danych klas (type: lista 3-elementowa)
        :return: np.array o rozmiarze p_table.shape[0] x len(lambda_we) -> 7x3
        """
        ilosc_kolejek = p_table.shape[1] - 2  # odejmujemy wejscie i wyjscie
        A = np.zeros((ilosc_kolejek, ilosc_kolejek))
        np.fill_diagonal(A, 1)
        b = lambda_we * p_table[0][1:-1]
        for row_i, p_row in enumerate(p_table):
            if row_i >= 1 and row_i <= A.shape[0]:
                A[row_i - 1] = A[row_i - 1] - p_table[row_i][1:-1]

        A = np.transpose(A)
        return np.linalg.solve(A, b)  # rozwiazywanie ukladu rownan

    @staticmethod
    def _single_Q_ir(i, r, m_i, ro_i, ro_ir):
        """
        :param i: nr systemu, type: int 0-6
        :param r: nr klasy, type: int   0-2
        :param m_i: liczba stanowisk w każdym systemie, type: list[int], len: 7
        :param ro_i: względna intensywność obsługi w i-tym systemie, ro_i <=1, type:  list[double], len: 7
        :param ro_ir: względna intensywność obsługi r-tej klasy w i-tym systemie, type: np.array
        :return: wartość Q_ir dla danego systemu i klasy
        """

        if np.isnan(ro_ir[i][r]):
            return 0
        else:
            x = ro_ir[i][r] / (1 - ro_i)
            y = (m_i[i] * ro_i) ** m_i[i] / (np.math.factorial(m_i[i]) * (1 - ro_i))

            temp_list = []
            for ki in range(m_i[i] - 1):
                if ki > 0:
                    temp_list.append((m_i[i] * ro_i) ** ki / np.math.factorial(ki))

            z = 1 / (sum(temp_list) + (m_i[i] * ro_i) ** m_i[i] / np.math.factorial(m_i[i]) * (1 / (1 - ro_i)))

            return x * z * y

    def f_cost_full(self, m_i):
        """
        Funkcja kosztów. W przypadku gdy nie są spełnione war ergodyczności zwraca 1000

        :param m_i: liczba stanowisk w każdym systemie, type: list[int], len: 7
        :return: wartość funkcji, type: double
        """
        m_i = [int(m_ii) for m_ii in m_i]
        r = 3
        lambdas_ir_list = []
        p_tables_r = [self.p_shortdist, self.p_longdist, self.p_transport]
        for p_table, lambda_we in zip(p_tables_r, self.lambda_we_r):
            # print("---------------")
            # print("lambda = ", lambda_we)
            lambdas_ir = self._calculate_lambdas_ir(p_table, lambda_we)

            # print("Przepustowość tej klasy lambda_ir = ", lambdas_ir)
            # print("\n")

            lambdas_ir_list.append(lambdas_ir)

        # lambdas_ir_arr

        lambdas_ir_arr = np.zeros((7, 3))
        for i in range(7):
            for r in range(3):
                lambdas_i = lambdas_ir_list[r]
                if r <= 1:
                    # pasażerskie
                    if i < 5:
                        lambdas_ir_arr[i][r] = lambdas_i[i]
                    elif i == 5:
                        lambdas_ir_arr[6][r] = lambdas_i[i]

                elif r == 2:
                    # towarowe
                    if i < 4:
                        lambdas_ir_arr[i][r] = lambdas_i[i]
                    elif i >= 4 and i < 6:
                        lambdas_ir_arr[i + 1][r] = lambdas_i[i]

        self.lambdas_ir_arr = lambdas_ir_arr
        # print("lambdas_ir_arr")
        # print(lambdas_ir_arr)


        # ro_iR
        ro_ir = np.zeros([7, 3])
        for i in range(7):
            for r in range(3):
                ro_ir[i][r] = lambdas_ir_arr[i][r] / (m_i[i] * self.mi_ir[i][r])

        self.ro_ir = ro_ir

        for i in range(7):
            ro_i = np.nan_to_num(ro_ir)[i].sum()
            if ro_i >= 1:
                # print("Warunek ro_i < 1 nie spełniony")
                return 1000

        # print("ro_ir:")
        # print(ro_ir)

        # liczba niezajetych kanałow
        m_nzi = []
        # sprawdzenie war. ergodyczności
        for i in range(7):
            ro_i = np.nan_to_num(ro_ir)[i].sum()
            m_nzi.append(m_i[i] - m_i[i] * ro_i)

        # print("m_nzi:")
        # print(m_nzi)

        self.m_nzi = m_nzi

        # Q_ir

        Q_ir = np.zeros((7, 3))

        for i in range(7):
            ro_i = np.nan_to_num(ro_ir)[i].sum()
            for r in range(3):
                Q_ir[i][r] = self._single_Q_ir(i, r, m_i, ro_i, ro_ir)

        self.Q_ir = Q_ir

        K_ir = np.zeros((7, 3))

        for i in range(7):
            for r in range(3):
                K_ir[i][r] = m_i[i] * ro_ir[i][r] + Q_ir[i][r]

        self.K_ir = K_ir

        # sumowanie kosztów
        costs = np.zeros((len(m_nzi), r))
        for i in range(len(m_nzi)):
            for j in range(r):
                costs[i][j] = self.C1_ij[i][j] * Q_ir[i][j] + self.C2_i[i] * m_nzi[i]

        return sum(sum(np.nan_to_num(costs)))

    def print_results(self, mis):
        cost_val = self.f_cost_full(mis)

        print(f"Rozwiązanie {mis}")

        print(f"Koszt: {cost_val}")

        print("\nprzepustowość r-tej klasy w i-tym systemie")
        print(self.lambdas_ir_arr)

        print("\n względna intensywność obsługi r-tej klasy w i-tym systemie")
        print(self.ro_ir)

        print("\n Q_ir średnia liczba zgłoszeń r-tej klasy w kolejce i-tego systemu")
        print(self.Q_ir)

        print("\n K_ir średnia liczba zgłoszeń r-tej klasy w i-tym systemie")
        print(self.K_ir)

        print("\n liczba niezajętych kanałów w systemie i-tym")
        print(self.m_nzi)





siec = AirportNetwork()

# m_i = [2,  # 0 pasy ladowania
#        6,  # 1 przeglad
#        2,  # 2 serwisowe stanowska
#        5,  # 3 tankowanie
#        3,  # 4 pasazerowie
#        2,  # 5 towary
#        3]  # 6 pasy wylotu
sample_m_i = [1, 5., 3, 2., 4., 1., 2.]

res = siec.f_cost_full(sample_m_i)
print(f"\n\n Result for m_i = {sample_m_i}: \n", res)
siec.print_results(sample_m_i)


ograniczenia = np.array([[1, 2],  # 0 pasy ladowania
                         [1, 5],  # 1 przeglad
                         [1, 5],  # 2 serwisowe stanowska
                         [1, 5],  # 3 tankowanie
                         [1, 5],  # 4 pasazerowie
                         [1, 2],  # 5 towary
                         [1, 3]])  # wylot
ga_model = ga(function=siec.f_cost_full, dimension=7, variable_type='int',
              variable_boundaries=ograniczenia,
              algorithm_parameters={'max_num_iteration': 100,
                                    'population_size': 10,
                                    'mutation_probability': 0.05,
                                    'elit_ratio': 0.2,
                                    'crossover_probability': 0.8,
                                    'parents_portion': 0.3,
                                    'crossover_type': 'uniform',
                                    'max_iteration_without_improv': None})

for it in range(5):


    ga_model.run()

    print("\n --------------------------")
    myList = ga_model.report
    minIndex = myList.index(min(myList))
    print(f"Rozwiązanie znalezione w {minIndex} iteracji")

    mi_opt = ga_model.best_variable

    siec.print_results(mi_opt)


    print("Kontynuować? ")
    input()