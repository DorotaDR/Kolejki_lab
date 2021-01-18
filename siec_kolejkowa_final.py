import numpy as np
from geneticalgorithm import geneticalgorithm as ga





class AirportNetwork:

    def __init__(self):

        lambda_we_r = [1.8, 1.2, 1]

        mi_ir = np.zeros([7, 3])
        mi_ir[0, :] = np.array([4, 4, 4])
        mi_ir[1, :] = np.array([1, 0.8, 1])
        mi_ir[2, :] = np.array([0.25, 0.2, 0.3])
        mi_ir[3, :] = np.array([3, 3 / 2, 2])
        mi_ir[4, :] = np.array([3 / 2, 1, 0])
        mi_ir[5, :] = np.array([0, 0, 2])
        mi_ir[6, :] = np.array([2, 2, 2])

        _columns = ["Wejście",
                    "Lądowanie",
                    "Przegląd",
                    "Serwis",
                    "Tankowanie",
                    "Pasażerowie",
                    "Wylot",
                    "Wyjscie"]

        p_shortdist = np.zeros((len(_columns), len(_columns)))
        p_shortdist[0][1] = 1
        p_shortdist[1][2] = 1
        p_shortdist[2][3] = 0.1
        p_shortdist[2][4] = 0.2
        p_shortdist[2][5] = 0.7
        p_shortdist[3][2] = 1
        p_shortdist[4][5] = 1
        p_shortdist[5][6] = 1
        p_shortdist[6][7] = 1

        p_longdist = np.zeros((len(_columns), len(_columns)))
        p_longdist[0][1] = 1
        p_longdist[1][2] = 1
        p_longdist[2][3] = 0.1
        p_longdist[2][4] = 0.9
        p_longdist[3][2] = 1
        p_longdist[4][5] = 1
        p_longdist[5][6] = 1
        p_longdist[6][7] = 1

        p_transport = np.zeros((len(_columns), len(_columns)))
        p_transport[0][1] = 1
        p_transport[1][2] = 1
        p_transport[2][3] = 0.1
        p_transport[2][4] = 0.7
        p_transport[2][5] = 0.2
        p_transport[3][2] = 1
        p_transport[4][5] = 1
        p_transport[5][6] = 1
        p_transport[6][7] = 1

        ## C1_ij
        # TODO set
        C1_ij = np.ones((7, 3))

        ## C2_i
        # TODO set
        C2_i = np.ones(7)

    #TODO: Inne kontruktory ustawianie kosztow C1, C2


    def _calculate_lambdas_ir(self, p_table, lambda_we):
        ilosc_kolejek = p_table.shape[1] - 2  # odejmujemy wejscie i wyjscie
        A = np.zeros((ilosc_kolejek, ilosc_kolejek))
        np.fill_diagonal(A, 1)
        b = lambda_we * p_table[0][1:-1]
        for row_i, p_row in enumerate(p_table):
            if row_i >= 1 and row_i <= A.shape[0]:
                A[row_i - 1] = A[row_i - 1] - p_table[row_i][1:-1]

        A = np.transpose(A)
        return np.linalg.solve(A, b)

    def _single_Q_ir(self, i, r, m_i, ro_i, ro_ir):
        x = ro_ir[i][r] / (1 - ro_i)
        y = (m_i[i] * ro_i) ** m_i[i] / (np.math.factorial(m_i[i]) * (1 - ro_i))

        temp_list = []
        for ki in range(m_i[i] - 1):
            if ki > 0:
                temp_list.append((m_i[i] * ro_i) ** ki / np.math.factorial(ki))

        z = 1 / (sum(temp_list) + (m_i[i] * ro_i) ** m_i[i] / np.math.factorial(m_i[i]) * (1 / (1 - ro_i)))

        return x * z * y






    def f_cost_full(self, m_i):
        m_i = [int(m_ii) for m_ii in m_i]
        r = 3
        lambdas_ir_list = []
        p_tables_r = [self.p_shortdist, self.p_longdist, self.p_transport]
        for p_table, lambda_we in zip(p_tables_r, self.lambda_we_r):
            print("---------------")
            print("lambda = ", lambda_we)
            lambdas_ir = self._calculate_lambdas_ir(p_table, lambda_we)

            print("Przepustowość tej klasy lambda_ir = ", lambdas_ir)
            print("\n")

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

        print("lambdas_ir_arr")
        print(lambdas_ir_arr)

        # ro_iR
        ro_ir = np.zeros([7, 3])
        for i in range(7):
            for r in range(3):
                ro_ir[i][r] = lambdas_ir_arr[i][r] / (m_i[i] * self.mi_ir[i][r])

        print("ro_ir:")
        print(ro_ir)


        ## liczba niezajetych kanałow
        m_nzi = []
        for i in range(7):
            ro_i = np.nan_to_num(ro_ir)[i].sum()
            if ro_i > 1:
                print("Warunek ro_i < 1 nie spełniony")
                return 100
            else:
                m_nzi.append(m_i[i] - m_i[i] * ro_i)

        print("m_nzi:")
        print(m_nzi)


        ## Q_ir
        Q_ir = np.zeros((7, 3))

        for i in range(7):
            ro_i = np.nan_to_num(ro_ir)[i].sum()
            for r in range(3):
                Q_ir[i][r] = self._single_Q_ir(i, r, m_i, ro_i, ro_ir)



        costs = np.zeros((len(m_nzi), r))
        for i in range(len(m_nzi)):
            for j in range(r):
                costs[i][j] = self.C1_ij[i][j] * Q_ir[i][j] + self.C2_i[i] * m_nzi[i]

        return sum(sum(np.nan_to_num(costs)))

siec = AirportNetwork()

# m_i = [2,  # 0 pasy ladowania
#        6,  # 1 przeglad
#        2,  # 2 serwisowe stanowska
#        5,  # 3 tankowanie
#        3,  # 4 pasazerowie
#        2,  # 5 towary
#        3]  # 6 pasy wylotu
m_i =  [1., 4., 1., 1, 2., 2., 1.]
# len(m_i)

res = siec.f_cost_full(m_i)
print("\n\n Restult: \n", res)

ograniczenia = np.array([[1, 2],  # 0 pasy ladowania
                         [1, 5],  # 1 przeglad
                         [1, 5],  # 2 serwisowe stanowska
                         [1, 5],  # 3 tankowanie
                         [1, 5],  # 4 pasazerowie
                         [1, 2],  # 5 towary
                         [1, 3]])  #
ga_model = ga(function=siec.f_cost_full, dimension=7, variable_type='int',
              variable_boundaries=ograniczenia,
              algorithm_parameters={'max_num_iteration': 100,
                                    'population_size': 5,
                                    'mutation_probability': 0.05,
                                    'elit_ratio': 0.2,
                                    'crossover_probability': 0.8,
                                    'parents_portion': 0.3,
                                    'crossover_type': 'uniform',
                                    'max_iteration_without_improv': None})

for it in range(5):
    ga_model.run()
