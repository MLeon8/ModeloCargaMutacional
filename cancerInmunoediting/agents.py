import mesa
import numpy as np
import random 

class CancerCell(mesa.Agent):
    def __init__(self, unique_id, model, mu, sigma,pos):
        super().__init__(unique_id, model)
        self.antiTumor = False
        self.age = 1
        self.pos = pos
        
    def step(self):
        # print("Activación CC")
        self.age += 1    

class InIScell(mesa.Agent):
    def __init__(self, unique_id, model, mu, sigma, maxAge,pos):
        super().__init__(unique_id, model)
        self.antiTumor = True
        self.maxAge = maxAge
        self.age = 1
        self.pos = pos

class CellM(InIScell):
    def __init__(self, unique_id, model, mu, sigma, mu2, sigma2, maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
        self.prProTumor = np.random.normal(mu/(mu+mu2),sigma)
        self.mu2 = mu2
        self.sigma2 = sigma2
        self.maxAge = {'M1':maxAge[0],'M2':maxAge[1]}
        self.CellType = ''
        
    def polarization(self):
        '''
        Solo se va a polarizar cuando existan células tumorales, no necesariamente
        cuando detecte antígeno
        '''
        if self.CellType == '':
            CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
            if len(CCs)>0 and np.random.uniform(0,100) < self.model.successAttackM:
                # print("Polarizando M")
                CC = self.random.choice(CCs)
                self.model.grid.place_agent(self,CC.pos)
                if np.random.uniform(0,100) < self.model.changeM1ToM2:
                    self.CellType = 'M1'
                    self.model._contM1Cells += 1
                else:
                    self.CellType = 'M2'
                    self.model._contM2Cells += 1
                self.model._contMCells -= 1
        
    def tumorInteraction(self):
        if (interac:=(np.random.uniform(0,100) < self.model.successAttackM1)) and self.CellType == 'M1':
            CC = [elem for elem in self.model.grid.get_cell_list_contents((16,16)) if elem.__class__.__name__ == 'CancerCell']
            randomCC = random.sample(CC, min(len(CC), self.model.MaxDeactivatingCCByM1))
            for cell in randomCC:
                cell.age += 1
        elif self.CellType == 'M2' and not(interac):
            CC = [elem for elem in self.model.grid.get_cell_list_contents((16,16)) if elem.__class__.__name__ == 'CancerCell']
            randomCC = random.sample(CC, min(len(CC), self.model.MaxActivatingCCByM1))
            for cell in randomCC:
                cell.age -= 1
    
    def die(self):
        probKill=4
        if self.CellType != "" and self.maxAge[self.CellType] < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            if self.CellType == "M1":
                self.model._contM1Cells -= 1
            else:
                self.model._contM2Cells -= 1
        self.age += 1
    
    def step(self):
        # print("Activación M")
        self.polarization()
        self.tumorInteraction()
        self.die() 
        
class CellN(InIScell):
    '''
    Los neutrófilos no tienen la capacidad para atacar a las células tumorales
    '''
    def __init__(self, unique_id, model, mu, sigma, mu2, sigma2, maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
        self.prProTumor = np.random.normal(mu/(mu+mu2),sigma)
        self.mu2 = mu2
        self.sigma2 = sigma2
        self.maxAge = {'N1':maxAge[0],'N2':maxAge[1]}
        self.CellType = ''
        
    def polarization(self):
        '''
        Solo se va a polarizar cuando existan células tumorales, no necesariamente
        cuando detecte antígeno
        '''
        if self.CellType == '':
            CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
            if len(CCs)>0 and np.random.uniform(0,100) > self.model.successAttackN:
                CC = self.random.choice(CCs)
                self.model.grid.place_agent(self,CC.pos)
                if np.random.uniform(0,100) < self.model.changeM1ToM2:
                    self.CellType = 'N1'
                    self.model._contN1Cells += 1
                else:
                    self.CellType = 'N2'
                    self.model._contN2Cells += 1
                self.model._contNCells -= 1
        
    def tumorInteraction(self):
        if (interac:=(np.random.uniform(0,100) < self.model.successAttackN1)) and self.CellType == 'N1':
            CC = [elem for elem in self.model.grid.get_cell_list_contents((16,16)) if elem.__class__.__name__ == 'CancerCell']
            randomCC = random.sample(CC, min(len(CC), self.model.MaxDeactivatingCCByM1))
            for cell in randomCC:
                cell.age += 0.1
        elif self.CellType == 'N2' and not(interac):
            CC = [elem for elem in self.model.grid.get_cell_list_contents((16,16)) if elem.__class__.__name__ == 'CancerCell']
            randomCC = random.sample(CC, min(len(CC), self.model.MaxActivatingCCByM1))
            for cell in randomCC:
                cell.age -= 0.1 
                
    def die(self):
        probKill=4
        if  self.CellType != "" and self.maxAge[self.CellType] < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            if self.CellType == "N1":
                self.model._contN1Cells -= 1
            else:
                self.model._contN2Cells -= 1
        self.age += 1
    
    def step(self):
        # print("Activación N")
        self.polarization()
        self.tumorInteraction()
        self.die() 

class CellNK(InIScell):
    def __init__(self, unique_id, model, mu, sigma,maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
    
    def attack(self):
        # Aumentar x2 
        killPercent = 5 
        # CancerCell = TypeVar("CancerCell")
        # print(f"CC clases : {isinstance(CancerCell, type)}")
        # CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
        # print(f"CCs: {len(CCs)}")
        if self.model._contCancerCells > 0:
            CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
            CC = self.random.choice(CCs)
            self.model.grid.place_agent(self,CC.pos)
        # CCs = [elem for elem in self.model.grid.get_cell_list_contents([self.pos]) if elem.__class__.__name__ == 'CancerCell']
        # print(f"NK células NK {len(CCs)}, su posición es {self.pos}")
            if np.random.uniform(0,100) < self.model.successAttackNk:
                # print("NK atacando")
                if self.age > (self.model.maxAgeM2 + self.model.maxAgeN2)/2 and np.random.uniform(0,100) < killPercent :
                    # CC = self.random.choice(CCs)
                    self.model.grid.remove_agent(CC)
                    self.model.schedule.remove(CC)
                    self.model._muertesCCporNK += 1
                    self.model._contCancerCells -= 1
                # CellN = TypeVar("CellN", bound=InIScell)
                # print(f"N2 clases : {isinstance(CellN,InIScell)}")
                # CellsN2 = [agent for agent in self.model.schedule.agents_by_type[CellN].values() if agent.CellType == 'N2' ]
                # print(f"N2: {CellsN2}")
                if self.model._contN2Cells > 0:
                    cellsN = [agent for agent in self.model.schedule.agents_by_type[CellN].values()]
                    cellN = self.random.choice(cellsN)
                    # print(f"N: {len(cellsN)}")
                    if cellN.CellType == "N2" and cellN.age > self.model.maxAgeN2 and np.random.uniform(0,100) < killPercent :
                        self.model.grid.remove_agent(cellN)
                        self.model.schedule.remove(cellN)
                        self.model._contN2Cells -= 1
                # CellM = TypeVar("CellM", bound=InIScell)
                # CellsM2 = [agent for agent in self.model.schedule.agents_by_type[CellM].values() if agent.CellType == 'M2' ]
                # print(f"M2: {CellsM2}")
                if self.model._contM2Cells > 0:
                    cellsM = [agent for agent in self.model.schedule.agents_by_type[CellM].values() ]
                    cellM = self.random.choice(cellsM)
                    # print(f"M: {len(cellsM)}")
                    if cellM.CellType == "M2" and cellM.age > self.model.maxAgeM2 and np.random.uniform(0,100) < killPercent :
                        self.model.grid.remove_agent(cellM)
                        self.model.schedule.remove(cellM)
                        self.model._contM2Cells -= 1
    
    def die(self):
        probKill=4
        if self.maxAge < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.model._contNKCells -= 1
        self.age += 1
    
    def step(self):
        # print("Activación NK")
        self.attack()
        self.die()

class CellGammaDeltaT(InIScell):
    def __init__(self, unique_id, model, mu, sigma,maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
    
    def attack(self):
        killPercent = 5 
        # # CancerCell = TypeVar("CancerCell")
        # # print(f"CC clases : {isinstance(CancerCell, type)}")
        # # CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
        # # print(f"CCs: {len(CCs)}")
        if self.model._contCancerCells > 0:
            CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
            CC = self.random.choice(CCs)
            self.model.grid.place_agent(self,CC.pos)
        # CCs = [elem for elem in self.model.grid.get_cell_list_contents([self.pos]) if elem.__class__.__name__ == 'CancerCell']
        # print(f"NK células NK {len(CCs)}, su posición es {self.pos}")
            if np.random.uniform(0,100) < self.model.successAttackGDT:
                # print("NK atacando")
                if self.age > (self.model.maxAgeM2 + self.model.maxAgeN2)/2 and np.random.uniform(0,100) < killPercent :
                    # CC = self.random.choice(CCs)
                    self.model.grid.remove_agent(CC)
                    self.model.schedule.remove(CC)
                    self.model._muertesCCporNK += 1
                    self.model._contCancerCells -= 1
    
    def die(self):
        probKill=4
        if self.maxAge < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.model._contNKCells -= 1
        self.age += 1
    
    def step(self):
        # print("Activación NK")
        self.attack()
        self.die()
        
class dendriticCells(InIScell):
    def __init__(self, unique_id, model, mu, sigma,maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
    
    def die(self):
        probKill=4
        if self.maxAge < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.model._contDCells -= 1
        self.age += 1
    
    def recruitAIS(self):
        if self.model.probNeoantigeno*random.randrange(0,100) > random.uniform(0,100):
            numCells = self.model._contTCells + self.model._contThCells
            
            newCells = self.model.recruitIS(self.model.schedule.time, numCells,self.model.aIS, self.model.isGrowthFactor)
            for i in range(newCells):
                pos = self.model.randomCoordinates()
                cell = self.model.TCell(self.model.next_id(), self.model, self.model.mu1, self.model.sigma1, self.model.maxAgeT,pos)
                self.model.schedule.add(cell)
                self.model.grid.place_agent(cell, pos)
            self._conTCells += newCells    
            # verificar si es necesario realizar otra función para reclutar, 
            # la cual divide entre 2.6 
            
            for i in range(newCells):
                pos = self.model.randomCoordinates()
                cell = self.model.CellTh(self.model.next_id(), self.model, self.model.mu1, self.model.sigma1, self.model.mu2, self.model.sigma2, self.model.maxAgeTh,pos)
                self.model.schedule.add(cell)
                self.model.grid.place_agent(cell, pos)
            self._contThCells += newCells   
        
class AdIScell(mesa.Agent):
    def __init__(self, unique_id, model, mu, sigma, maxAge,pos):
        super().__init__(unique_id, model)
        self.antiTumor = True
        self.age = 1
        self.maxAge = maxAge
        self.pos = pos

class TCell(AdIScell):
    def __init__(self, unique_id, model, mu, sigma,maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
        self.successAttackT = self.model.successAttackT
    
    def die(self):
        probKill=4
        if self.maxAge < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.model._contTCells -= 1
        else:
            self.age += 1
    
    def Attack(self):
        '''
        En la implementación de NetLogo se permite que esta célula mate a más de 
        una célula por vez, esto debido a que ejecuta este método para todas las células
        tumorales que se encuentren en la misma posición
        Sin embargo, la cantidad de células 
        '''
        # Aumentar x2 
        killPercent = 4
        
        CCs = [ agent for agent in self.model.schedule.agents_by_type[CancerCell].values() ]
        if len(CCs) and np.random.uniform(0,100) < self.model.successAttackT:
            CC = self.random.choice(CCs)    
            self.model.grid.place_agent(self,CC.pos)
        CCs = [elem for elem in self.model.grid.get_cell_list_contents(self.pos) if elem.__class__.__name__ == 'CancerCell']
        for agent in CCs:
            if np.random.uniform(0,100) < self.successAttackT and self.age > (self.model.maxAgeM2 + self.model.maxAgeN2)/2 and np.random.uniform(0,100) < killPercent:
                self.model.grid.remove_agent(agent)
                self.model.schedule.remove(agent)
                self.model._contCancerCells -= 1
                self.model._muertesCCporT += 1
                
    def step(self):
        # print("Activación T")
        self.Attack()
        self.die()
    
class ThCell(AdIScell):
    def __init__(self, unique_id, model, mu, sigma,maxAge,pos):
        super().__init__(unique_id, model, mu, sigma,maxAge,pos)
        self.CellType = ""
        
    def polarization(self):
        if self.CellType == '':
            if prob := np.random.uniform(0,100) < self.model.changeTh[1] :
                self.CellType = "Th1"
                self.model._contTh1 += 1
            elif prob < self.model.changeTh[2]:
                self.CellType = "Th2"
                self.model._contTh2 += 1
            elif prob < self.model.changeTh[3]:
                self.CellType = "Th17"
                self.model._contTh17 += 1
            else:
                self.CellType = "Treg"
                self.model._contTreg += 1
            self.model._contTCells -= 1
    
    def interaction(self):
        if self.CellType == "Th1":
            self.model.prRecuit = max(self.model.prRecuit+1,100)
            self.model.successAttackT = max(self.model.successAttackT+1, 100)
        elif self.CellType == "Th2":
            self.model.changeT[1] = min(self.model.changeT[2]-5,max(self.model.changeT[0]+5,self.model.changeT[1]-1))
            # Preguntar cómo es que las Th2 afectan a la polarización de los macrófagos
        elif self.CellType == "Th17":
            self.model.changeN1ToN2 = max(self.model.changeN1ToN2+1,100)
        elif self.CellType == "Threg":
            pass
        # aquí te quedaste
        
    def strengthening(self):
        TCs = [ agent for agent in self.model.schedule.agents_by_type[TCell].values() ]
        if len(TCs) and np.random.uniform(1,100) < self.model.successAttackThToT:
            TC = self.random.choice(TCs)
            self.model.grid.place_agent(self,TC.pos)
        TCs = [elem for elem in self.model.grid.get_cell_list_contents(self.pos) if elem.__class__.__name__ == 'TCell']
        for TC in TCs:
            TC.successAttackT += 0.1
        
    def die(self):
        probKill=4
        if self.maxAge < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.model._contThCells -= 1
        else:
            self.age += 1
    
    def step(self):
        # print("Activación Th")
        self.strengthening()
        self.die()

class TregCell(AdIScell):
    '''
    En el programa de NetLogo una célula Treg puede fortalecer a más de una célula 
    T o Th, esto se debe a que se considera la parte geométrica
    Estas células modifican la edad de Th y T, pero como estas no dependen de la edad,
    no tienen un impacto dentro del programa
    '''
    def __init__(self, unique_id, model, mu, sigma, maxAge,pos):
        super().__init__(unique_id, model, mu, sigma, maxAge,pos)
            
    def strengthening(self):
        if np.random.uniform(0,100) < 50:
            TCs = [ agent for agent in self.model.schedule.agents_by_type[TCell].values() ]
            if len(TCs) and np.random.uniform(1,100) < self.model.successAttackTregToT:
                TC = self.random.choice(TCs)
                self.model.grid.place_agent(self,TC.pos)
                TC.age -= 0.1
        else:
            ThCs = [ agent for agent in self.model.schedule.agents_by_type[ThCell].values() ]
            if len(ThCs) and np.random.uniform(1,100) < self.model.successAttackTregToTh:
                ThC = self.random.choice(ThCs)
                self.model.grid.place_agent(self,ThC.pos)
        TCs = [elem for elem in self.model.grid.get_cell_list_contents(self.pos) if elem.__class__.__name__ == 'TCell']
        for TC in TCs:
            TC.age -= 0.1
        ThCs = [elem for elem in self.model.grid.get_cell_list_contents(self.pos) if elem.__class__.__name__ == 'ThCell']
        for ThC in ThCs:
            ThC.age += 0.1
    
    def die(self):
        probKill=4
        if self.maxAge < self.age and np.random.uniform(0,100) < probKill:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.model._contTregCells -= 1
        else:
            self.age += 1
    
    def step(self):
        # print("Activación Treg")
        self.strengthening()
        self.die()