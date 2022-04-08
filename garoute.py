
import numpy
from random import *
import matplotlib.pyplot as plt

class Ambulance:
    def __init__(self, x=None, y=None, PNo=None ):

        self.x = x
        self.y = y
        self.plateNo = PNo

    def getX(self):
        return self.x

    def getY(self):
        return self.y


    def plateNo(self):
        return self.plateNo


class Emergency :

    def __init__(self, x=None, y=None, ENo=None):
        self.x = x
        self.y = y
        self.emergencyNo = ENo

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def emergencyNo(self):
        return self.emergencyNo

    def distance(self, location):
        xDis = abs(self.x - location.x)
        yDis = abs(self.y - location.y)
        distance = numpy.sqrt((xDis ** 2) + (yDis ** 2))
        return distance

class Hospital :

    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def distance(self, location):
        xDis = abs(self.x - location.x)
        yDis = abs(self.y - location.y)
        distance = numpy.sqrt((xDis ** 2) + (yDis ** 2))
        return distance


class Location:

    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def distance(self, location):
        xDis = abs(self.x - location.x)
        yDis = abs(self.y - location.y)
        distance = numpy.sqrt((xDis ** 2) + (yDis ** 2))
        return distance

    def distance(self, location):
        xDis = abs(self.x - location.x)
        yDis = abs(self.y - location.y)
        distance = numpy.sqrt((xDis ** 2) + (yDis ** 2))
        return distance



class LocationManager:

    locations = []
    def addLocation(self, location):
        self.locations.append(location)

    def getLocation(self, index):
        return self.locations[index]

    def numberOfLocations(self):
        return len(self.locations)

class HospitalManager:

    hospitals = []

    def addHospital(self, hospital):
        self.hospitals.append(hospital)

    def getHospital(self, index):
        return self.hospitals[index]

    def numberOfHospitals(self):
        return len(self.hospitals)

class EmergencyManager:

    emergencies = []

    def addEmergency(self, emergency):
        self.emergencies.append(emergency)

    def getEmergency(self, index):
        return self.emergencies[index]

    def numberOfEmergencies(self):
        return len(self.emergencies)

class Chromosome:

    def __init__(self,emergency, locationmanager, hospitalmanager  ):
        self.hospitalManager = hospitalmanager
        self.locationManager = locationmanager
        self.emergency = emergency
        self.fitness = 0.0
        self.distance = 0
        self.locations = []
        self.LocationGeneCount = 0



    def generateIndividual(self):

        numOfLocations = self.locationManager.numberOfLocations()

        randomGeneSize = randint(2, numOfLocations-1  )

        for i in range(0, randomGeneSize):
         self.locations.append(None)
        for i in range (0, randomGeneSize):
            self.locations[i] = self.locationManager.getLocation(i)

        shuffle(self.locations)

        self.locations[0] = self.emergency

        for i in range(randomGeneSize, numOfLocations -1 ):
         self.locations.append(None)

        self.locations.append(None)
        randomIndex = randint(0,self.hospitalManager.numberOfHospitals() - 1)
        randomHospital = self.hospitalManager.getHospital( randomIndex )
        self.locations[numOfLocations -1] = randomHospital


        self.fitness = 0.0
        self.distance = 0
        if(self.LocationGeneCount==0):
            self.LocationGeneCount = randomGeneSize -1

    def getChromosomeSize(self):

        return len(self.locations)


    def getDistance(self):
        if self.distance == 0:
            locationsDistance = 0

            locationIndex =0
            while(self.locations[locationIndex] is not None and locationIndex < (len(self.locations)-1)):

                if(self.locations[locationIndex+1] is not None):
                    fromLocation = self.locations[locationIndex]
                    destinationLocation = self.locations[locationIndex+1]
                    locationsDistance += fromLocation.distance(destinationLocation)

                locationIndex = locationIndex + 1

            fromLocation = self.locations[locationIndex-1]
            destinationLocation = self.locations[-1]
            locationsDistance += fromLocation.distance(destinationLocation)


            self.distance = locationsDistance
            return self.distance

        else:
            return self.distance

    def getDistance2(self):
        if self.distance == 0:
            locationsDistance = 1
            if(self.LocationGeneCount==1):
                locationsDistance = self.locations[0].distance(self.locations[1])

            else:
                locationIndex =0
                while(self.locations[locationIndex] is not None and locationIndex <(len(self.locations)-1)):

                    if(self.locations[locationIndex+1] is not None):
                        fromLocation = self.locations[locationIndex]
                        destinationLocation = self.locations[locationIndex+1]
                        locationsDistance += fromLocation.distance(destinationLocation)

                    locationIndex = locationIndex + 1

                # fromLocation = self.locations[locationIndex-1]
                # destinationLocation = self.locations[-1]
                # locationsDistance += fromLocation.distance(destinationLocation)


            self.distance = locationsDistance
            return self.distance

        else:
            return self.distance


    def getFitness(self):

        if self.fitness == 0:
            f = float(self.locations[0].distance(self.locations[-1]))
            f1 = float(self.getDistance2())

            self.fitness =1/(f*f1)
            return self.fitness
        else:
            return self.fitness

    def DuplicatePresent(self,index,location):

        for i in range (1,index+1):
            if(self.locations[i] == location):
                return True

        else:
            return False


class Population:
    def __init__(self, emergency, locationmanager, hospitalmanager, populationSize, initialise):
        self.chromosomes = []
        for i in range(0, populationSize):
            self.chromosomes.append(None)

        if initialise:
            for i in range(0, populationSize):
                newchromosome = Chromosome(emergency, locationmanager, hospitalmanager )
                newchromosome.generateIndividual()
                self.saveChromosome(i, newchromosome)

    def __setitem__(self, key, value):
        self.chromosomes[key] = value

    def __getitem__(self, index):
        return self.chromosomes[index]

    def saveChromosome(self, index, chromosome):
        self.chromosomes[index] = chromosome


    def getFittest(self):
        fittest = self.chromosomes[0]
        for i in range(1, self.populationSize()):
            f1 = fittest.getFitness()
            f2 = self.__getitem__(i).getFitness()
            if f1 <= f2:
                fittest = self.chromosomes[i]
        return fittest

    def getLowestFittestIndex(self):
        fittest = self.chromosomes[0]
        index = 0
        for i in range(1, self.populationSize()):
            if fittest.getFitness() >= self.__getitem__(i).getFitness():
                fittest = self.chromosomes[i]
                index = i
        return index


    def populationSize(self):
        return len(self.chromosomes)





class RoutePlannerGA:

    def __init__(self, locationmanager, hospitalmanager, emergency):
        self.locationManager = locationmanager
        self.hospitalManager = hospitalmanager
        self.Emergency = emergency
        self.mutationRate = 0.7
        self.tournamentSize = 10
        self.elitism = True
        self.cc = Chromosome(self.Emergency, self.locationManager,self.hospitalManager)
        self.cc2 = Chromosome(self.Emergency, self.locationManager,self.hospitalManager)
        self.locationPos1 = 0
    def evolvePopulation(self, pop):
        newPopulation = Population(self.Emergency, self.locationManager , self.hospitalManager, pop.populationSize(), False)
        elitismOffset = 0
        if self.elitism:
            newPopulation.saveChromosome(0, pop.getFittest())
            elitismOffset = 1

        for i in range(elitismOffset, newPopulation.populationSize()):
            parent1 = self.tournamentSelection(pop)
            parent2 = self.tournamentSelection(pop)
            if (parent1 == parent2):
                while (parent1 == parent2):
                    parent2 = self.tournamentSelection(pop)
            child = self.crossover(parent1, parent2)
            newPopulation.saveChromosome(i, child)

        # for i in range(elitismOffset, newPopulation.populationSize()):
        #     self.mutate(newPopulation.chromosomes[i])
        return newPopulation

    def evolvePopulation2(self, pop):
        index= pop.getLowestFittestIndex()
        parent1 = self.tournamentSelection(pop)
        parent2 = self.tournamentSelection(pop)

        if(parent1 == parent2 ):
            while(parent1==parent2):
                parent2 = self.tournamentSelection(pop)
        child1 = self.crossover(parent1, parent2)
        self.cc = child1
        child2 = self.crossover(parent2, parent1)
        self.cc2 = child2
        self.mutate(child1)

        self.mutate(child2)
        child1.getFitness()
        child2.getFitness()
        pop.chromosomes[index] = child1
        index = pop.getLowestFittestIndex()


    def crossover(self, parent1, parent2):
        child = Chromosome(self.Emergency, self.locationManager, self.hospitalManager)

        child.locations.append(child.emergency)
        index = 1

        # for i in range(0, self.locationManager.numberOfLocations()):
        #     child.locations.append(None)

        half1 = int(parent1.LocationGeneCount/2)
        if(half1 < 1):
            half1 = 1

        parent1RandIndex = int(randint(1,half1) )#from start to index into child of parent1's
        half2 = int(parent2.LocationGeneCount/2)
        if (half2 < 1):
            half2 = 1
        parent2RandIndex = int(randint(half2  , parent2.LocationGeneCount ) ) ##from start to index into child of parent2's

        childGeneCount = (parent1RandIndex + ((parent2.LocationGeneCount+1) - parent2RandIndex ))
        minusCount = 0
        if(childGeneCount <= (parent1.getChromosomeSize() - 1)): #to check for the variable length gene count

            for i in range (1,parent1RandIndex+1):
                child.locations.append(parent1.locations[i]) #append upto first parents index

            for i in range (parent1RandIndex+1, childGeneCount+1):
                if(child.DuplicatePresent(parent1RandIndex,parent2.locations[parent2RandIndex])): #repairs for duplicate value
                    minusCount +=1
                else:
                    child.locations.append(parent2.locations[parent2RandIndex])  # append upto first parents index
                    parent2RandIndex+=1
            childGeneCount -= minusCount
            for i in range (childGeneCount+1, parent1.getChromosomeSize()-1):
                child.locations.append(None)

            randHospital = random() # if 1 selects parent1's hospital
            if(randHospital > 0.1):
                child.locations.append(parent1.locations[-1])
            else:
                child.locations.append(parent2.locations[-1])

        child.LocationGeneCount = childGeneCount
        return child

    def mutate(self, chromosome):
        if chromosome.LocationGeneCount == 1:
           return

        else:
            for self.locationPos1 in range(1, chromosome.LocationGeneCount):
                r = random()
                if r < self.mutationRate:
                    locationPos2 = randint(1, chromosome.LocationGeneCount )

                    loc1 = chromosome.locations[self.locationPos1]
                    loc2 = chromosome.locations[locationPos2]

                    chromosome.locations[locationPos2] = loc1
                    chromosome.locations[self.locationPos1] =loc2






    def tournamentSelection(self, pop):

        tournament = Population(self.Emergency,self.locationManager, self.hospitalManager,self.tournamentSize, False)

        for i in range(0, self.tournamentSize):
            randomId = int(random() * pop.populationSize())

            tournament.saveChromosome(i, pop.__getitem__(randomId))

        fittest = tournament.getFittest()

        return fittest



if __name__ == '__main__':

    locationManager = LocationManager()
    hospitalManager = HospitalManager()
    emergency = Emergency(7,6)

    chromosomeSize = 35
    NoOfHos = 5
    PopulationSize = 50

    for i in range(0,chromosomeSize):
        location = Location(randint(0,60), randint(0,60))
        locationManager.addLocation(location)

    hospital = Hospital(0,0)
    hospitalManager.addHospital(hospital)

    hospital = Hospital(30, 30)
    hospitalManager.addHospital(hospital)
    hospital = Hospital(5, 10)
    hospitalManager.addHospital(hospital)
    hospital = Hospital(0, 60)
    hospitalManager.addHospital(hospital)
    hospital = Hospital(50, 50)
    hospitalManager.addHospital(hospital)

    # Initialize population
    pop = Population(emergency, locationManager, hospitalManager, 100, True )
    print ("Initial bestFitness: " + str(pop.getFittest().getFitness()))

    # initialize genetic algorithm for route planning
    ga = RoutePlannerGA(locationManager, hospitalManager, emergency)
    # pop = ga.evolvePopulation(pop)
    ga.evolvePopulation2(pop)

    for i in range(0, 1000):
         # pop = ga.evolvePopulation(pop)
         ga.evolvePopulation2(pop)


        # Print final results
    fittest = pop.getFittest()
    print("Finished")
    print("Final distance: " + str(pop.getFittest().getDistance()))
    print ("Solution:")
    print (pop.getFittest().getFitness())



    #--------------------------GRAPH PLOTTING------------------------------------
    x = []
    y = []
    xloc = []
    yloc = []
    index = 0
    while(pop.getFittest().locations[index] is not None):
        x.append(pop.getFittest().locations[index].getX())
        y.append(pop.getFittest().locations[index].getY())
        index+=1
    x.append(pop.getFittest().locations[-1].getX())
    y.append(pop.getFittest().locations[-1].getY())

    for i in range (0 , chromosomeSize):
        xloc.append(locationManager.locations[i].getX())
        yloc.append(locationManager.locations[i].getY())

    plt.scatter(xloc,yloc,edgecolors='darkblue')
    plt.scatter(emergency.getX(), emergency.getY(), edgecolors='red')

    xhos = []
    yhos = []
    for i in range(0, NoOfHos):
        xhos.append(hospitalManager.hospitals[i].getX())
        yhos.append(hospitalManager.hospitals[i].getY())

    plt.scatter(xhos, yhos, edgecolors='green')

    plt.plot(x, y, color='lightblue', linewidth=3)


    plt.show()
