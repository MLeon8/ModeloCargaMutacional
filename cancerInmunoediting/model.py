import mesa
import numpy as np
import random
from cancerInmunoediting.scheduler import RandomActivationByTypeFiltered
from cancerInmunoediting.agents import CancerCell, CellNK, CellM, CellN,\
    TCell, ThCell, TregCell, dendriticCells, CellGammaDeltaT, dendriticCells
'''
NOTAS:
    TODO lo que esté muestreado de una función de probabilidad gaussiana debe de 
    aplicarsele el abs()
'''
class CancerInmunoediting(mesa.Model):
    description = (
        "A model for simulating cancer inmunoediting, quantified by Ising Hamiltonian"
        )
    verbose = False
    
    
    def __init__(self,meanIS, stdIS, meanCancer, stdCancer, ml=1000):
        super().__init__()
        self.sigma1 = stdIS
        self.mu1 = meanIS
        self.sigma2 = stdCancer
        self.mu2 = meanCancer
        self.HAntiCancer = 0
        self.HProCancer = 0
        self.HTME = 0
        
        self.nProCancer = 0
        self.nAntiCancer = 0
        
        self.wMatrix = np.zeros([8,8],dtype=float)
        self.qMatrix = np.zeros([8,8],dtype=float)

        # rango de 0.001 a  1000
        # ¿se debe de normalizar (dividir entre 1000)?
        # ml asociado al cáncer
        # para un tipo de cáncer específico
        
        # el usuario ingresa el tipo de cáncer
        # self.type_cancer = type_cancer
        self.type_cancer_mu = 1
        self.type_cancer_sigma = 100
        

        # La proporción debe de variar dentro de un rango dado por el usuario
        self.probNeoantigeno = self.cargaMutacional*random.randrange(60,90)*0.01 
        
        self.qMatrix = [
            [0,1,1,1,1,1,1,1],
            [1,0,1,1,0,0,0,0],
            [1,1,0,1,0,0,0,0],
            [1,1,1,0,0,0,0,0],
            [1,0,0,0,0,1,1,0],
            [1,0,0,0,1,0,1,0],
            [1,0,0,0,1,1,0,0],
            [1,0,0,0,0,0,0,0]
            ]
        
        self.xVector = np.zeros(8,dtype=float)
        
        self._contCancerCells = 0
        self._contNKCells = 0
        self._contNCells = 0
        self._contMCells = 0
        self._contM1Cells = 0
        self._contM2Cells = 0
        self._contN1Cells = 0
        self._contN2Cells = 0
        self._contDCells = 0
        self._contTCells = 0
        self._contTregCells = 0
        self._contThCells = 0
        self._contTh1 = 0
        self._contTh2 = 0
        self._contTh17 = 0
        self._contTreg = 0
        
        self.GrowthFactor =  0
        
        self.width = 32
        self.height = 32
        
        self.schedule = RandomActivationByTypeFiltered(self)
        self.grid = mesa.space.MultiGrid(self.width , self.height, True)
        # print(
        #     [
        #         self.schedule.time,
        #         self.schedule.get_type_count(CancerCell),
        #         self._contCancerCells,
        #         self.schedule.get_type_count(CellNK),
        #         self._contNKCells,
        #         self.schedule.get_type_count(CellM, lambda x: x.CellType == ''),
        #         self._contMCells,
        #         self.schedule.get_type_count(CellN, lambda x: x.CellType == ''),
        #         self._contNCells,
        #         self.schedule.get_type_count(CellM, lambda x: x.CellType == 'M1'),
        #         self._contM1Cells,
        #         self.schedule.get_type_count(CellN, lambda x: x.CellType == 'N1'),
        #         self._contN1Cells,
        #         # self.schedule.get_type_count(TCell),
        #         # self.schedule.get_type_count(ThCell),
        #         # self.schedule.get_type_count(TregCell),
        #         self.schedule.get_type_count(CellM, lambda x: x.CellType == 'M2'),
        #         self._contM2Cells,
        #         self.schedule.get_type_count(CellN, lambda x: x.CellType == 'N2'),
        #         self._contN2Cells,
        #         # self._contNKCells + self._contNCells - self._contN2Cells \
        #         #     + self._contMCells - self._contM2Cells,
        #         self._contNKCells + self._contN1Cells + self._contM1Cells,
        #         self._contCancerCells + self._contM2Cells + self._contN2Cells,
        #         # self.HAntiCancer,
        #         # self.HProCancer,
        #         # self.HTME
        #     ])
        self.dictFuction = {'CC' : lambda uniqueID, model, mu, sigma,k ,t0 : CancerCell(uniqueID, model, mu, sigma, k, t0)}
        self.datacollector = mesa.DataCollector(
            {
                # "CancerCells" : lambda m : m.schedule.get_type_count(CancerCell),
                "CancerCells" : lambda m : m._contCancerCells,
                # "CellsNK": lambda m: m.schedule.get_type_count(CellNK),
                "NK": lambda m: m._contNKCells,
                # "CellsM1": lambda m: m.schedule.get_type_count(CellM, lambda x: x.CellType == 'M1'),
                "M1": lambda m: m._contM1Cells,
                # "CellsM2": lambda m: m.schedule.get_type_count(CellM, lambda x: x.CellType == 'M2'),
                "M2": lambda m: m._contM2Cells,
                # "CellsN1": lambda m: m.schedule.get_type_count(CellN, lambda x: x.CellType == 'N1'),
                "N1": lambda m: m._contN1Cells,
                # "CellsN2": lambda m: m.schedule.get_type_count(CellN, lambda x: x.CellType == 'N2'),
                "N2": lambda m: m._contN2Cells,
                "T": lambda m: m._contTCells,
                "Th": lambda m: m._contThCells,
                "Treg": lambda m: m._contTregCells,
                # "AntiCancer": lambda m : m.schedule.get_count( lambda x: x.antiTumor),
                # Esta expresión es distinta a la reportada, debido a que de otra
                # forma la población es negativa
                "AntiCancer": lambda m: m._contNKCells + m._contNCells + m._contN1Cells \
                    + m._contMCells + m._contM1Cells, 
                # "ProCancer" : lambda m: m.schedule.get_count( lambda x: not(x.antiTumor)),
                "ProCancer": lambda m: m._contCancerCells + m._contM2Cells + m._contN2Cells,
                "HAntiCancer" : lambda m: m.HAntiCancer,
                "HProCancer" : lambda m: m.HProCancer,
                "HTME": lambda m: m.HTME,
                "TumorGrowthRate": lambda m: m.GrowthFactor
                }
            )
        
        self.initialCancerCells = 0
        self.initialNaturalKillers = 0
        self.initialMacrofagues = 0
        self.initialNeutrophils = 0
        self.initialTCells = 0
        self.initialThCells = 0
        self.initialTregCells = 0
        self.initialDendriticCells = 0
        
        self.maxAge = 10
        self.maxCells = 100
        self.maxCellsRecruit = 100
        self.maxCellsActivate = 100
        self.maxCellsDeactivate = 100
        
        # solo debe de haber 2 ataques, los de NK y T
        self.successAttackN = 0
        self.successAttackM = 0
        self.successAttackT = 0
        self.successAttackTregToT = 0
        self.successAttackTregToTh = 0
        self.successAttackThToT = 0
        self.successAttackNk = 0
        self.successAttackTh = 0
        self.successAttackM1 = 0
        self.successAttackN1 = 0
        self.successAttackM2 = 0
        self.successAttackN2 = 0
        
        self.contadorMuerteIS = 0
        self.contadorMuerteCC = 0
        
        self.maxAgeM1 = 0
        self.maxAgeN1 = 0
        self.maxAgeM2 = 0
        self.maxAgeN2 = 0
        self.maxAgeNK = 0
        self.maxAgeT = 0
        self.maxAgeTh = 0
        self.maxAgeTreg = 0
        self.maxAgeDC = 0
        
        self.changeM1ToM2 =  0
        self.changeN1ToN2 = 0
        
        self.changeThToTh1 = 25
        self.changeThToTh2 = 25
        self.changeThToTh17 = 25
        self.changeThToTreg = 25
        
        self.changeTh = [0,25,50,75,100]
        
        
        self.MaxDeactivatingCCByM1 = 0
        self.MaxDeactivatingCCByN1 = 0

        self.MaxActivatingCCByM1 = 0
        self.MaxActivatingCCByN1 = 0
        
        self.isGrowthFactor = 0
        self.aIS = 0
        
        self.ccGrowthFactor = 0
        self.aCC = 0
        
        self.running = True
        self.datacollector.collect(self)
        self.updateXVector()
        self.initializeData()      
        self.initWMatrix()
        
        for i in range(self.initialCancerCells):
            # x   y
            # -16 16 
            # x : 32 + xNETLOGO
            # y : yNETLOGO
            pos = (16,16)
            cell = CancerCell(self.next_id(), self, self.mu2, self.sigma2, pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
        
        for i in range(self.initialNaturalKillers):
            pos = self.randomCoordinates()
            cell = CellNK(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeNK,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
            
        for i in range(self.initialMacrofagues):
            pos = self.randomCoordinates()
            cell = CellM(self.next_id(), self, self.mu1, self.sigma1, self.mu2, self.sigma2, [self.maxAgeM1, self.maxAgeM2],pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
        
        for i in range(self.initialNeutrophils):
            pos = self.randomCoordinates()
            cell = CellN(self.next_id(), self, self.mu1, self.sigma1, self.mu2, self.sigma2, [self.maxAgeN1, self.maxAgeN2],pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
            
        for i in range(self.initialTCells):
            pos = self.randomCoordinates()
            cell = TCell(self.next_id(), self, self.mu1, self.sigma1,self.maxAgeT,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
            
        for i in range(self.initialThCells):
            pos = self.randomCoordinates()
            cell = ThCell(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeTh,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
        
        for i in range(self.initialTregCells):
            pos = self.randomCoordinates()
            cell = TregCell(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeTreg,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
            
        for i in range(self.initialDendriticCells):
            pos = self.randomCoordinates()
            cell = dendriticCells(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeDC,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
        
        
        self._muertesCCporNK = 0
    
    def randomCoordinates(self): # report cordinates[sign crd]
        x = random.randint(0,self.width-1)
        y = random.randint(0,self.height-1)
        return (x,y)
    
    def initializeData(self):
        self.initialCancerCells = abs(int(np.random.normal(self.mu2,self.sigma2)*100))
        self.initialNaturalKillers = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.initialMacrofagues = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.initialNeutrophils = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.initialTCells = int(np.random.normal(self.mu1,self.sigma1)*100)
        self.initialThCells = int(np.random.normal(self.mu1,self.sigma1)*100)
        self.initialTregCells = int(np.random.normal(self.mu1,self.sigma1)*100)
        self.initialDendriticCells = int(np.random.normal(self.mu1,self.sigma1)*100)
        
        self._contCancerCells = abs(self.initialCancerCells)
        self._contNKCells = abs(self.initialNaturalKillers)
        self._contMCells = abs(self.initialMacrofagues)
        self._contNCells = abs(self.initialNeutrophils)
        self._contDCells = abs(self.initialDendriticCells)
        self._contTCells = abs(self.initialTCells)
        self._contThCells = abs(self.initialThCells)
        self._contTregCells = abs(self.initialTregCells)
        
        self.successAttackN = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackM = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackT = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackTregToT = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackTregToTh = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackThToT = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackNk = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackTh = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackM1 = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackN1 = abs(int(np.random.normal(self.mu1,self.sigma1)*100))
        self.successAttackM2 = abs(int(np.random.normal(self.mu2,self.sigma2)*100))
        self.successAttackN2 = abs(int(np.random.normal(self.mu2,self.sigma2)*100))
        
        self.maxAgeM1 = abs(int(np.random.normal(self.mu1,self.sigma1)*10))
        self.maxAgeN1 = abs(int(np.random.normal(self.mu1,self.sigma1)*10))
        self.maxAgeM2 = abs(int(np.random.normal(self.mu2,self.sigma2)*10)) 
        self.maxAgeN2 = abs(int(np.random.normal(self.mu2,self.sigma2)*10)) 
        self.maxAgeNK = abs(int(np.random.normal(self.mu1,self.sigma1)*10))
        self.maxAgeT = abs(int(np.random.normal(self.mu1,self.sigma1)*10))
        self.maxAgeTh = abs(int(np.random.normal(self.mu1,self.sigma1)*10))
        self.maxAgeTreg = abs(int(np.random.normal(self.mu1,self.sigma1)*10))
        self.maxAgeDC = abs(int(np.random.normal(self.mu1,self.sigma1)*10))

        self.changeM1ToM2 =  abs(int(np.random.normal(self.mu1/(self.mu1 + self.mu2),self.sigma1)*100))
        self.changeN1ToN2 = abs(int(np.random.normal(self.mu1/(self.mu1 + self.mu2),self.sigma1)*100))
        
        self.MaxDeactivatingCCByM1 = abs(int(self.maxCellsDeactivate * np.random.normal(self.mu1,self.sigma1)))
        self.MaxDeactivatingCCByN1 = abs(int(self.maxCellsDeactivate * np.random.normal(self.mu1,self.sigma1)))

        self.MaxActivatingCCByM1 = abs(int(self.maxCellsActivate * np.random.normal(self.mu2,self.sigma2)))
        self.MaxActivatingCCByN1 = abs(int(self.maxCellsActivate * np.random.normal(self.mu2,self.sigma2)))
        
        self.isGrowthFactor = abs(np.random.normal(self.mu1,self.sigma1)*100)
        self.aIS = self.derLogistic(10*(self.isGrowthFactor * 0.01),50)
        
        self.ccGrowthFactor = abs(np.random.normal(self.mu2,self.sigma2)*100)
        self.aCC = self.derLogistic(10*(self.ccGrowthFactor * 0.01),50)
        
        # iniciales = [
        #     4,
        #      71,
        #      68,
        #      37,
        #      66,
        #      59,
        #      41,
        #      #4,
        #      #71,
        #      #68,
        #      #37,
        #      97,
        #      0,
        #      152,
        #      96,
        #      24,
        #      59,
        #      11,
        #      25,
        #      72,
        #      67,
        #      4,
        #      26,
        #      4,
        #      3,
        #      10,
        #      4,
        #      7,
        #      7,
        #      2,
        #      5,
        #      179,
        #      122,
        #      128,
        #      12,
        #      77,
        #      38,
        #      9.635328315748348,
        #      0.013732375032412758,
        #      61.95419835480327,
        #      0.03799592737219082
        #     ]
        
        # self.initialCancerCells = iniciales[0]
        # print(self.initialCancerCells)
        # self.initialNaturalKillers = iniciales[1]
        # self.initialMacrofagues = iniciales[2]
        # self.initialNeutrophils = iniciales[3]
        # self.initialTCells = iniciales[4]
        # self.initialThCells = iniciales[5]
        # self.initialTregCells = iniciales[6]
        
        # self._contCancerCells = self.initialCancerCells
        # self._contNKCells = self.initialNaturalKillers 
        # self._contMCells = self.initialMacrofagues
        # self._contNCells = self.initialNeutrophils 
        
        # self.successAttackN = iniciales[7]
        # self.successAttackM = iniciales[8]
        # self.successAttackT = iniciales[9]
        # self.successAttackTregToT = iniciales[10]
        # self.successAttackTregToTh = iniciales[11]
        # self.successAttackThToT = iniciales[12]
        # self.successAttackNk = iniciales[13]
        # self.successAttackTh = iniciales[14]
        # self.successAttackM1 = iniciales[15]
        # self.successAttackN1 = iniciales[16]
        # self.successAttackM2 = iniciales[17]
        # self.successAttackN2 = iniciales[18]
        
        # self.maxAgeM1 = iniciales[19]
        # self.maxAgeN1 = iniciales[20]
        # self.maxAgeM2 = iniciales[21]
        # self.maxAgeN2 = iniciales[22]
        # self.maxAgeNK = iniciales[23]
        # self.maxAgeT = iniciales[24]
        # self.maxAgeTh = iniciales[25]
        # self.maxAgeTreg = iniciales[26]
        
        # self.changeM1ToM2 =  iniciales[27]
        # self.changeN1ToN2 = iniciales[28]
        
        # self.MaxDeactivatingCCByM1 = iniciales[29]
        # self.MaxDeactivatingCCByN1 = iniciales[30]

        # self.MaxActivatingCCByM1 = iniciales[31]
        # self.MaxActivatingCCByN1 = iniciales[32]
        
        # self.isGrowthFactor = iniciales[33]
        # self.aIS = iniciales[34]
        
        # self.ccGrowthFactor = iniciales[35]
        # self.aCC = iniciales[36]
        
        # print(f""" 
        #       ------------------
        #       {self.initialCancerCells},
        #       {self.initialNaturalKillers},
        #       {self.initialMacrofagues},
        #       {self.initialNeutrophils},
        #       {self.initialTCells},
        #       {self.initialThCells},
        #       {self.initialTregCells},
        #      # {self._contCancerCells},
        #      # {self._contNKCells},
        #      # {self._contMCells},
        #      # {self._contNCells},
        #       {self.successAttackN},
        #       {self.successAttackM},
        #       {self.successAttackT},
        #       {self.successAttackTregToT},
        #       {self.successAttackTregToTh},
        #       {self.successAttackThToT},
        #       {self.successAttackNk},
        #       {self.successAttackTh},
        #       {self.successAttackM1},
        #       {self.successAttackN1},
        #       {self.successAttackM2},
        #       {self.successAttackN2},
        #       {self.maxAgeM1},
        #       {self.maxAgeN1},
        #       {self.maxAgeM2},
        #       {self.maxAgeN2},
        #       {self.maxAgeNK},
        #       {self.maxAgeT},
        #       {self.maxAgeTh},
        #       {self.maxAgeTreg},
        #       {self.changeM1ToM2},
        #       {self.changeN1ToN2},
        #       {self.MaxDeactivatingCCByM1},
        #       {self.MaxDeactivatingCCByN1},
        #       {self.MaxActivatingCCByM1},
        #       {self.MaxActivatingCCByN1},
        #       {self.isGrowthFactor},
        #       {self.aIS},
        #       {self.ccGrowthFactor},
        #       {self.aCC}
        #       ---------------------
        #       """)
        
        
    
    def updateXVector(self):
        self.xVector = [
            -1 * self._contCancerCells, 
             # 1 * self._contNCells - self._contN1Cells - self._contN2Cells,
             1 * self._contNCells,
             1 * self._contN1Cells,
            -1 * self._contN2Cells,
             # 1 * self._contMCells - self._contM1Cells - self._contM2Cells,
             1 * self._contMCells ,
             1 * self._contM1Cells,
            -1 * self._contM2Cells,
             1 * self._contNKCells
            ]
        
    def initWMatrix(self):
        p12 = self.successAttackN * 0.01
        p23 = self.changeN1ToN2 * 0.01
        p24 = 1 - p23
        p13 = self.successAttackN1 * 0.01
        p14 = 1 - p13
        p15 = self.successAttackM * 0.01
        p56 = self.changeM1ToM2 * 0.01
        p57 = 1 - p56
        p16 = self.successAttackM1 * 0.01
        p17 = 1 -p16
        p18 = self.successAttackNk * 0.01
        self.wMatrix = [
            [0.0,p12,p13,p14,p15,p16,p17,p18],
            [p12,0.0,p23,p24,0.0,0.0,0.0,0.0],
            [p13,p23,0.0,0.5,0.0,0.0,0.0,0.0],
            [p14,p24,0.5,0.0,0.0,0.0,0.0,0.0],
            [p15,0.0,0.0,0.0,0.0,p56,p57,0.0],
            [p16,0.0,0.0,0.0,p56,0.0,0.0,0.0],
            [p17,0.0,0.0,0.0,p57,0.0,0.0,0.0],
            [p18,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            ]
        
        
    # prototipo de funcion para generar celulas (No terminada)
    def createInitialCells(self,**args):
        for i in range(args['noCells']):
            cell = self.dictFuction['typeCell'](elem for elem in args)
            self.schedule.add(cell)
        
            
    def logistic(self,x,a,b):
        k=4
        res = []
        for i in x:
            res.append(1 + a/(1 + b*(np.e**(-k*(i-0.5)))))
        return res
    
    def derLogistic(self,a,b):
        x = np.arange(0,50)
        tickSpread = 50
        xx=x/(tickSpread-1)
        x1=self.logistic(xx,a,1)[1:]
        x2=self.logistic(xx,a,1)[0:-1]
        res = []
        for i in range(len(x1)):
            # print(x1[i],x2[i])
            try:
                res.append(x1[i]/x2[i])
            except:
                print(x1[i],x2[i])
        return max(res) - 1
    
    def gaussAprox(self,x,a):
        c=0.19 # desviación estándar 
        b=0.5  # valor máximo
        aux = a * ( np.e** ( - (( x-b )**2)/(2*(c**2))))
        # print(f"aux: {aux}, a: {a}")
        return aux#a*(np.e**(-(x-b)**2))/(2*c**2)
    
    def recruitCC(self, tiempo, celulasIniciales, a, growthFactor):
        # print("reclutamiento CC")
        tickSpread = 50
        dx = 1/(tickSpread -1 )
        recluta=self.gaussAprox(dx*tiempo, a)
        self.GrowthFactor = recluta
        # print(f"growthFactor :{recluta}, tiempo: {tiempo}")
        recluta = int(recluta)
        # print(f"recluta: {recluta}, celinicial:{celulasIniciales}, growth " )
        aux = growthFactor
        
        # print(f"aux: {aux}, x:{dx*tiempo}, a:{a}")
        try:
            result= celulasIniciales*recluta
            # print(f"recluta: {recluta}, result: {result}, celulasIniciales: {celulasIniciales}, GrowthFactor: {aux}")
            aux = np.random.uniform(0,100)
            
            if result == 0 and aux < growthFactor:
                result=3
            # print(f"recluta: {recluta}, result: {result}")
            # print(f"growthFactor: {growthFactor} \t muestra: {aux} \t Nuevas Células: {result}")
            # celulas.append(result+celulas[-1])
            
            return result
        except:
            print("error")
        
    def recruitIS(self, tiempo, celulasIniciales, a, growthFactor):
        # print("reclutamiento IS")
        tickSpread = 50
        dx = 1/(tickSpread -1 )
        recluta=int(self.gaussAprox(dx*tiempo, a))
        # print(f"recluta: {recluta}, celinicial:{celulasIniciales}, growth:{growthFactor} ")
        aux = growthFactor
        try:
            result= celulasIniciales*recluta//2.6
            # print(f"recluta: {recluta}, result: {result}, celulasIniciales: {celulasIniciales}, GrowthFactor: {aux}")
            aux= np.random.uniform(0,100)
            if result == 0 and aux < growthFactor:
                result=1
            # print(f"growthFactor: {growthFactor} \t muestra: {aux} \t Nuevas Células: {result}")

            # print(f"recluta: {recluta}, result: {result}")
            # celulas.append(result+celulas[-1])
            return int(result)
        except:
            print("error")
        
    def cancerGrowth(self):
        self.newCells = self.recruitCC(self.schedule.time, (self.schedule.get_type_count(CancerCell)),self.aCC, self.ccGrowthFactor)
        # print(f"Células máximas {self.aCC}")
        # print(f"Nuevas células {self.newCells}")
        self._contCancerCells += self.newCells
        
        for i in range(self.newCells):
            pos = (16,16)
            cell = CancerCell(self.next_id(), self, self.mu2, self.sigma2,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
        # print(self.newCells)
        return self.newCells
    
    def ISrecruit(self):        
        if self._contCancerCells:
            celulasIS = self.schedule.get_type_count(CellNK) + self.schedule.get_type_count(CellN) \
                + self.schedule.get_type_count(CellM, lambda x: x.antiTumor)
            
            newCells = self.recruitIS(self.schedule.time, celulasIS,self.aIS, self.isGrowthFactor)
            for i in range(int(newCells*0.5)):
                pos = self.randomCoordinates()
                cell = CellNK(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeNK,pos)
                self.schedule.add(cell)
                self.grid.place_agent(cell, pos)
            self._contNKCells += newCells    
            
            newCells = self.recruitIS(self.schedule.time, celulasIS,self.aIS, self.isGrowthFactor)
            for i in range(newCells):
                pos = self.randomCoordinates()
                cell = CellN(self.next_id(), self, self.mu1, self.sigma1, self.mu2, self.sigma2, [self.maxAgeN1, self.maxAgeN2],pos)
                self.schedule.add(cell)
                self.grid.place_agent(cell, pos)
            self._contNCells += newCells   
            
            newCells = self.recruitIS(self.schedule.time, celulasIS,self.aIS, self.isGrowthFactor)
            for i in range(newCells):
                pos = self.randomCoordinates()
                cell = CellM(self.next_id(), self, self.mu1, self.sigma1, self.mu2, self.sigma2, [self.maxAgeM1, self.maxAgeM2],pos)
                self.schedule.add(cell)
                self.grid.place_agent(cell, pos)
            self._contMCells += newCells
        
    def AISrecuit(self):
        # rango de 0.001 a  1000
        # ¿se debe de normalizar (dividir entre 1000)?
        self.cargaMutacional = 1000
        # La proporción debe de variar dentro de un rango dado por el usuario
        self.probNeoantigeno = self.cargaMutacional*random.randrange(0,100)*.01
        self.probReclutamiendoDC = self.probNeoantigeno*random.uniform(60,90)*0.01
        
        for i in range(int(self.probReclutamiendoDC)):
            pos = self.randomCoordinates()
            cell = dendriticCells(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeNK,pos)
            self.schedule.add(cell)
            self.grid.place_agent(cell, pos)
        self._contNKCells += self.probReclutamiendoDC
        # probar para seno y próstata
        
        # if self.probReclutamiendoDC > random.uniform(0,100):
        #     celulasIS = self.schedule.get_type_count(CellNK) + self.schedule.get_type_count(CellN) \
        #         + self.schedule.get_type_count(CellM, lambda x: x.antiTumor)
            
        #     newCells = self.recruitIS(self.schedule.time, celulasIS,self.aIS, self.isGrowthFactor)
        #     for i in range(newCells):
        #         pos = self.randomCoordinates()
        #         cell = CellNK(self.next_id(), self, self.mu1, self.sigma1, self.maxAgeNK,pos)
        #         self.schedule.add(cell)
        #         self.grid.place_agent(cell, pos)
        #     self._contNKCells += newCells    
            
        #     newCells = self.recruitIS(self.schedule.time, celulasIS,self.aIS, self.isGrowthFactor)
        #     for i in range(newCells):
        #         pos = self.randomCoordinates()
        #         cell = CellN(self.next_id(), self, self.mu1, self.sigma1, self.mu2, self.sigma2, [self.maxAgeN1, self.maxAgeN2],pos)
        #         self.schedule.add(cell)
        #         self.grid.place_agent(cell, pos)
        #     self._contNCells += newCells   
            
        #     newCells = self.recruitIS(self.schedule.time, celulasIS,self.aIS, self.isGrowthFactor)
        #     for i in range(newCells):
        #         pos = self.randomCoordinates()
        #         cell = CellM(self.next_id(), self, self.mu1, self.sigma1, self.mu2, self.sigma2, [self.maxAgeM1, self.maxAgeM2],pos)
        #         self.schedule.add(cell)
        #         self.grid.place_agent(cell, pos)
        #     self._contMCells += newCells
        
    def updateH(self):
        # print(f"W: {self.wMatrix}")
        self.updateXVector()
        wAux = np.zeros([8,8])
        for i in range(len(self.wMatrix)):
            for j in range(len(self.wMatrix[0])):
                wAux[i][j] = self.wMatrix[i][j] * self.qMatrix[i][j]
        # print(f"WAux: {wAux}")
        # print(f"self.wMatrix: {self.wMatrix}")
        # print(f"self.qMatrix: {self.qMatrix}")
        self.HAntiCancer = 0
        self.HProCancer = 0
        self.HTME = 0
        aux = 0
        for i in range(len(self.xVector)):
            for j in range(len(self.xVector)):
                if abs(self.xVector[i]) > abs(self.xVector[j]):
                    aux = wAux[i][j] * self.xVector[i] * abs(self.xVector[j])
                else:
                    aux = wAux[i][j] * self.xVector[j] * abs(self.xVector[i])
                if aux >= 0:
                    self.HAntiCancer -= 0.5 * aux
                else:
                    self.HProCancer -= 0.5 * aux
        # print(f"HAnti: {self.HAntiCancer}, HPro: {self.HProCancer }")
        self.HTME = self.HAntiCancer + self.HProCancer

    def step(self):
        self._muertesCCporNK = 0
        self._muertesCCporT = 0
        self.ISrecruit()
        self.AISrecuit()
        self.updateH()
        # self.schedule.step()
        # CancerCell, CellNK, CellM, CellN, TCell, ThCell, TregCell
        self.schedule.step_type(CancerCell, True) 
        self.cancerGrowth()
        self.schedule.step_type(CellN, True) 
        
        self.schedule.step_type(CellM, True) 
        self.schedule.step_type(CellNK, True) 
        self.schedule.step_type(TregCell, True)
        self.schedule.step_type(ThCell, True)
        self.schedule.step_type(TCell, True)
        self.schedule.time += 1
        self.schedule.steps += 1
        
        self.datacollector.collect(self)
        if self.verbose:
            print(
                f"""
                    T  : {self.schedule.time},
                    CC : {self.schedule.get_type_count(CancerCell)},
                    {self._contCancerCells},
                    NK : {self.schedule.get_type_count(CellNK)},
                    {self._contNKCells},
                    M  : {self.schedule.get_type_count(CellM, lambda x: x.CellType == '')},
                    {self._contMCells},
                    N  : {self.schedule.get_type_count(CellN, lambda x: x.CellType == '')},
                    {self._contNCells},
                    M1 : {self.schedule.get_type_count(CellM, lambda x: x.CellType == 'M1')},
                    {self._contM1Cells},
                    N1 : {self.schedule.get_type_count(CellN, lambda x: x.CellType == 'N1')},
                    {self._contN1Cells},
                    M2 : {self.schedule.get_type_count(CellM, lambda x: x.CellType == 'M2')},
                    {self._contM2Cells},
                    N2 : {self.schedule.get_type_count(CellN, lambda x: x.CellType == 'N2')},
                    {self._contN2Cells},
                    Pro: {self._contNKCells + self._contNCells - self._contN2Cells + self._contMCells - self._contM2Cells},
                    {self._contNKCells + self._contN1Cells + self._contM1Cells},
                    Anti:{self._contCancerCells + self._contM2Cells + self._contN2Cells}
                """
            )
            # ProCancer": lambda m: m._contCancerCells + m._contM2Cells + m._contN2Cells,
            # print(
            #     [
            #         {self.schedule.time},
            #         {self.schedule.get_type_count(CancerCell)},
            #         {self._contCancerCells},
            #         {self.schedule.get_type_count(CellNK)},
            #         {self._contNKCells},
            #         {self.schedule.get_type_count(CellM, lambda x: x.CellType == '')},
            #         {self._contMCells},
            #         {self.schedule.get_type_count(CellN, lambda x: x.CellType == '')},
            #         {self._contNCells},
            #         {self.schedule.get_type_count(CellM, lambda x: x.CellType == 'M1')},
            #         {self._contM1Cells},
            #         {self.schedule.get_type_count(CellN, lambda x: x.CellType == 'N1')},
            #         {self._contN1Cells},
            #         # self.schedule.get_type_count(TCell),
            #         # self.schedule.get_type_count(ThCell),
            #         # self.schedule.get_type_count(TregCell),
            #         {self.schedule.get_type_count(CellM, lambda x: x.CellType == 'M2')},
            #         {self._contM2Cells},
            #         {self.schedule.get_type_count(CellN, lambda x: x.CellType == 'N2')},
            #         {self._contN2Cells},
            #         {self._contNKCells + self._contNCells - self._contN2Cells \
            #             + self._contMCells - self._contM2Cells},
            #         {self._contNKCells + self._contN1Cells + self._contM1Cells},
            #         {self._contCancerCells + self._contM2Cells + self._contN2Cells},
            #         # self.HAntiCancer,
            #         # self.HProCancer,
            #         # self.HTME
            #     ]
            # )

    def run_model(self, step_count=200):
        # if self.verbose:
        #     print("Initial number CellM: ", self.schedule.get_type_count(CellM))
        #     print("Initial number CellM1: ", self.schedule.get_type_count(CellM, lambda x: x.antiTumor))
        #     print(
        #         "Initial number CellM2: ",
        #         self.schedule.get_type_count(CellM, lambda x: not(x.antiTumor)),
        #     )
        for i in range(step_count):
            self.step()
        # if self.verbose:
        #     print("")
        #     print("Final number CellM: ", self.schedule.get_type_count(CellM))
        #     print("Final number CellM1: ", self.schedule.get_type_count(CellM, lambda x: x.antiTumor))
        #     print(
        #         "Final number CellM2: ",
        #         self.schedule.get_type_count(CellM, lambda x: not(x.antiTumor)),
        #     )