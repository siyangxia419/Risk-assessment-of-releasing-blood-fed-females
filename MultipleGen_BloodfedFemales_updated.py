import numpy as np
import numexpr as ne
# import math as m
# import sys as sys
# import matplotlib.pyplot as plt
import time
start_time = time.time()

np.set_printoptions(threshold=np.inf)

N_mosquito = 1000  # Number of wild mosquito populations
N_release = N_mosquito * 0.1  # Number of released mosquitoes
N_simulation = 10  # Number of simulations in each run

modelnumber = 0  # use this to run on laptop
# modelnumber=int(sys.argv[1])          # for cluster

# An matrix with all parameter values: each row is a parameter set and each column is one parameter
parameterarray = np.loadtxt("SensitivityAnalysisBloodfed", delimiter=",")
# Choose the running parameter values for a specific run
currentparameters = parameterarray[modelnumber, :]

with np.printoptions(precision=3, suppress=True):
    print(currentparameters)  # Print the current paramter values

hi, FF, vcw, vcr, ovi, bitephase, EIP, mubites, releaseratio, a, b, s, c, sdovi, sdhost, sdbite, sdEIP, sdVCw, sdVCr, fed, name, VC1, VC2, VC3, VC4, VC5, VC6, VC7, VC8, VC9, VC10 = currentparameters
# Assign parameter values to variables for easy use in the rest of the script

# hi           = % of Humans Infectious (dafault = 0.05)
# FF           = age of the mosquitoes when released (default = 5, i.e. the released mosquitoes are five days old)
# vcw          = mean vector competence of the wild mosquitoes (default = 0.524)
# vcr          = mean vector competence of the released mosquitoes (default = 0.2)
# ovi          = mean duration of the oviposition period (default = 4 days)
# bitephare    = mean duration of the host-seeking period (default = 4 days)
# EIP          = mean duration of the extrinsic incubation period (the time between a mosquito acquires the virus and it becomes infectious, default = 11 days)
# mubites      = mean number of bites in each host-seeking period (default = 2)
# releaseratio = size of the released population relative to the wild population (default = 1, which represents releasing 10% of the wild population size)
# a, b, s, c   = parameters that define the mosquito age-dependent survival function
# sdovi        = standard deviation of the oviposition period (default = 1 days)
# sdhost       = standard deviation of the host-seeking period (default = 1 days)
# sdbite       = standard deviation of the number of bites in each host-seeking period (default = 1)
# sdEIP        = standard deviation of EIP (default = 1 days)
# sdVCw        = standard deviation of the wild mosquito VC (default = 0.125)
# sdVCr        = standard deviation of the released mosquito VC (default = 0.125)
# fed          = proportion of the released mosquitoes being fed before the release (default = 1)
#    we are proposing to release female mosquitoes after blood-feeding them in the lab, but practically some individuals might fail to feed in the lab
# name         = index of runs (for tracking different runs on the cluster, not part of the simulation but is used for writing simulation output)
# VC1 - VC10   = the VC of the wild mosquito population throughout the 10 generations, expressed as proportion to the initial wild population VC (i.e. the first element is always 1)

# Note: vector competence represents the mosquito's ability to transmit diseases. We assume each mosquitoes have a VC range between 0 and 1.

VC_decline = np.array([VC1, VC2, VC3, VC4, VC5, VC6, VC7, VC8, VC9, VC10])

# The following four lines are the change of VC across the 10 generations for four difference scenarios
# VC_decline=np.array ([1,0.959,0.921,0.886,0.853,0.822,0.794,0.767,0.743,0.719])       #ONLY MALES
# For only males, there're effectively 0 released mosq bc no males bite (only VC change)
# VC_decline=np.array ([1,0.928,0.865,0.810,0.761,0.718,0.680,0.646,0.615,0.588])        #Blood-fed Fem + Males
# VC_decline=np.array ([1,0.956,0.916,0.879,0.845,0.813,0.783,0.756,0.730,0.707])       #Unfed Fem + Males
# VC_decline=np.array ([1,1,1,1,1,1,1,1,1,1])       #baseline


# Function to calulate the mortality of mosquitoes depending on their age
def age_dependence(x, a1, b1, s1, c1, delta_t):  # hazard rate calculator
    # defines the hazard rate function
    # adjust c to change agerage life expectancy, c=0.022 is 21 days
    """
    EDIT: Local variables are more efficient
    a1 = 1.799999999999999951e-03
    b1 = 1.416000000000000036e-01
    s1 = 1.072999999999999954e+00
    c1 = 2.199999999999999872e-02
    """
    return ne.evaluate('exp(-((a1 * exp(b1 * x)) / (1 + (a1 * s1 / b1) * (exp(b1 * x) - 1)) + c1) * delta_t)')


age_wild = age_dependence(np.arange(0, 100, 0.1), a, b, s, c, 0.1)
age_wild[-1] = 0  # force the mosquito to die at the last time breakpoint
survival_wild = np.cumprod(age_wild)

age_release = age_wild[int(10*FF):]
survival_release = np.cumprod(age_release)


# Function to sample the life span of a mosquito
def age_of_death(F, survival_p):  # age of death calculator
    # F = age at release (0 for wild, 5 for lab)
    """
    Old Code:
    age_m = np.zeros([1000])  # 1/10 of a day as the time resolution, simulate for 100 days
    age_matrix = np.asarray(age_m)
    i = 0
    while i < 999:
        age_matrix[i] = age_dependence(i / 10)
        i += 1
        # calculates the age dependent mortality at each 10th of a day
    """
    age_m = np.arange(0, 100, 0.1)
    age_matrix = age_dependence(age_m)

    fref = int(F * 10)

    """
    Old code:
    for x in range(0, 999):
        if fref > x:
            age_matrix[x] = 1  # released mosquitoes have no mortality before being released
    """
    age_matrix[0:fref] = 1

    # reshape the vector to a 100x10 matrix, with each row represents a single day and each column represent the 10th of each day
    age_matrixx = np.reshape(age_matrix, (100, 10))

    # The next two lines calculate the probability of the mosquito surviving to certain age (with 0.1 days as time resolution)
    achance_survive = np.cumprod(age_matrixx, axis=0)
    chance_survive = np.reshape(achance_survive, (1000, 1))

    rando = np.random.random_sample(1)
    itemindex = np.where(chance_survive < rando)
    # calls only where the array of chance-probabilities < random number value
    aitemindex = np.asarray(itemindex)
    deathcell = aitemindex[0, 0]
    DEATH = (
        deathcell - fref) / 10  # the time of death (for the release mosquitoes, this is the time of death after being released)

    return DEATH


# Function to simulate one single mosquito
# We first simulate a sequence of 10 gonotrophic cycle and them use the time of death to clip the sequence to preserve only the parts when the mosquito is alive
def one_mosquito(hi, binary, vcmu, vcsigma):  # for binary, 1=release 0=wild

    mu, sigma = vcmu, vcsigma  # Vector Competence
    vc = np.random.normal(mu, sigma,
                          1)  # each individual mosquito has a constant VC randomly sampled from a normal distribution

    mu, sigma = ovi, sdovi  # Duration of Oviposition
    do = np.random.normal(mu, sigma,
                          10)  # 10 numbers to represent 10 gonotrophic cycles (each includes a host-seeking and an oviposition period)
    trunc_do = np.clip(do, 3,
                       30)  # assume that each oviposition period must be at least 3 days and no more than 30 days
    ado = np.asarray([trunc_do])

    # For each gonotrophic cycle, we start with the oviposition phase and then the host-seeking phase
    # For wild mosquitoes, it requires a two-day maturation before the first host-seeking, which equilvalent to a two-day oviposition period
    # For consistency with the release mosquitoes, which enters the simulation in the oviposition period, we consider this "two-day" period in wild mosquitoes as its first oviposition peirod
    # So the next few lines change the first oviposition peirod of the wild mosquitoes
    firstO = np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1])  # Activated only for the wild pop
    add_2 = np.array([2, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    ado0 = ado * firstO
    ado2 = np.asarray(ado0 + add_2)
    if binary == 1:

        final_ato = ado  # "ado" for normal model, "ado0" for non-bloodfed mosq <-------------key input

        # Released population
    else:
        final_ato = ado2
        # wild population
    atdo = final_ato.transpose()
    attdo = atdo.transpose()

    mu, sigma = bitephase, sdhost  # Duration of Bloodfeeding
    db = np.random.normal(mu, sigma, 10)
    # At least spending one day in blood-feeding and no more than 30 days
    trunc_db = np.clip(db, 1, 30)
    adb = np.asarray([trunc_db])
    atdb = np.transpose(adb)

    mu, sigma = mubites, sdbite  # Bites per gonotrophic cycle (i.e. each host-seeking period)
    bites = np.random.normal(mu, sigma, 10)
    nb = np.array([bites])
    anb = np.round(nb)  # rounded bc: bites must be integers
    trunc_anb = np.clip(anb, 1, 10)  # at least one bites per cycle and at most 10 bites
    atnb = trunc_anb.transpose()

    mu, sigma = EIP, sdEIP  # Extrinsic Incubation Period
    eip_draw = np.random.normal(mu, sigma, 1)  # Only need one value for each individual
    eip = np.clip(eip_draw, 5, 30)  # EIP at least 5 days and at most 30 days

    infection_status = np.zeros(
        [10])  # a vector to record whether the mosquito got infected during each gonotrophic cycle
    # cycle through the two gonotrophic cycle:
    for x in range(0, 9):
        bitesincyc = int(atnb[x])  # extract the number of bites in this gonotrophic cycle
        # unclear why int() is necessary here, other situations, arrays work
        randn = np.random.random_sample((bitesincyc))
        for i in np.nditer(randn):  # determine whether each bite infect the mosquito
            if i < (hi * vc):  # hi * vc is the probability of the mosquito got affected biting any human
                infection_status[x] = 1
                # tests random number against probability
        if np.any(infection_status):  # The mosquito only needs to be infected once
            break
            # for all 1s, replace this line w 'infection_status[x+1] = 1'
    IS = np.asarray([infection_status])
    tIS = IS.transpose()

    empty_list = np.zeros(
        10)  # create a counter vector to record the days since the mosquitoes got affected for each gonotrophic cycle: if counter > EIP, the mosquito becomes infectious
    empty_array = np.array([empty_list])
    counter = empty_array.transpose()
    for x in range(0, 9):  # Calculates initial counter input
        biteday = (tIS[x])
        # calls the infectious binary signal
        bitephaseday = (atdb[x])
        # calls the length of biting cycle when mosquito was infected
        randn = np.random.random_sample(())
        if biteday == 1:
            counter[
                x] = randn * bitephaseday  # Within the host-seeking perid in which the mosquito got infected, what is the time between the infectious bites and the end of the biting phase
            # (random)*(days) is used as a 'randbetween' function

    gono = atdb + atdo  # combining blood-feeding duration and oviposition duration for full Gono-cycle length
    agono = np.asarray(gono)

    # Count the number of days since the mosquito becomes infected for all gonotrophic cycles
    for x in range(0, 9):  # Counter Mechanism
        countb4 = counter[x]
        cycle_days = agono[x + 1]
        # for x=0, countb4 also equals 0
        for i in np.nditer(countb4):
            if i > 0:
                counter[x + 1] = countb4 + cycle_days

    EmptyStatus = np.zeros(10)
    empty_arrayy = np.array([EmptyStatus])
    EIPtest = empty_arrayy.transpose()

    for x in range(0, 10):  # Infectious calculator
        # the number of days since the mosquito becomes infected at gonotrophic cycle x
        full_counter = counter[x]
        Incubation = eip
        # EIP is randomly sampled above
        for i in np.nditer(full_counter):
            if i >= Incubation:
                # when counter > EIP, the mosquito becomes infectious (able to infect another human host)
                EIPtest[x] = 1

    risky_bites = np.array(EIPtest * atnb)
    # binary test to eliminate non-infectious bites

    # Number of humans infected by this mosquito in each gonotrophic cycle
    empty_humans = np.zeros(10)
    tempty_humans = empty_humans.transpose()
    infected_humans = np.asarray(tempty_humans)
    for x in range(0, 10):
        bitesincyc = int(risky_bites[x])  # number of bites after the mosquito became infectious
        mosq_roulette = np.random.binomial(bitesincyc, (1 - hi))
        # binomial test to make sure bites on pre-infected humans aren't counted
        for i in np.nditer(bitesincyc):
            infected_humans[x] = mosq_roulette

    # The next few lines redefine the breakpoint of gonotrophic cycles
    # original gonotrophic cycle starts with oviposition -> biting
    # the new cycle: biting -> oviposition
    celltoinsert = np.zeros(1)  # Mosquito biting age calculation
    o = np.asarray(celltoinsert)
    # add a zero-day biting phase which pairs with the first oviposition phase
    oo = np.concatenate((o, trunc_db))
    # shifts the bite cycle phase down to add to next ovi phase
    shifted_bitephase = np.asarray(oo)
    clipped_end = np.delete(shifted_bitephase, 10)
    # eliminates exta value at the end so that it can add with ovi phase
    agemidcycle = clipped_end + attdo
    # Mosquito age at end of each gonotrophic cycle (i.e. the end of oviposition phase)
    biting_age = np.cumsum(agemidcycle)
    abiting_age = np.asarray(biting_age)

    agetodie = np.zeros(10)
    if binary == 1:  # Links the dif DOD functions to binary, binary = 1 indicates that we are simulating a release mosquito
        F = FF
    else:
        F = 0
    DOOD = age_of_death(F)
    # selects the death age function, assumes 5 days old for the released

    for x in range(0, 10):  # eliminates all transmissions after death
        biteage = abiting_age[x]  # age of the mosquito at the end of a gonotrophic cycle
        if biteage <= DOOD:  # If the mosquito has not died
            agetodie[x] = 1  # makes a binary check 1=alive
            if x < 9:
                # ensures that x+1 isn't out of range
                # the mosquito enters the next gonotrophic cycle, which means it has the chance to bite more people
                agetodie[x + 1] = 1
    # After this loop, agetodie should indicate whether the mosquito is still alive at the beginning of next gonotrophic cycle

    scaretest = 0  # Result of whether mosquito was infectious before died
    EIPttest = EIPtest.transpose()
    InfAndAlive = EIPttest * agetodie
    if np.any(InfAndAlive) > 0:
        scaretest = 1

    # number of humans the mosquito infected (after it died, the agetodie = 0, which makes the mosquito infected no more humans)
    living_infections = agetodie * infected_humans
    aliving_infections = np.asarray(living_infections)  # output array of actual transmissions

    # mic_drop = np.sum(aliving_infections)

    # print("\n")
    # print('The Mosquitos Have Left the Cage.\n')
    # print("Vector Competence: {}\n".format(vc))
    # print("Duration of Oviposition:\n {} \n".format(attdo))
    # print("Duration of Bloodfeeding:\n {} \n".format(atdb))
    # print("Number of Bites Per Cycle:\n {} \n".format(atnb))
    # print("Mosquito Infected:\n {} \n".format(tIS))
    # print("Extrinsic Incubation Period:\n {} \n".format(eip))
    # print("Incubation Counter:\n {} \n".format(counter))
    # print("INFECTED MOSQUITO ON THE LOOSE!:\n {} \n".format(EIPtest))
    # print("HUMAN INFECTIONS:\n {} \n".format(infected_humans))
    # print("Mosquito Age at Start of Bloodfeeding Phase:\n {} \n".format(biting_age))
    # print("Death Binary:\n {} \n".format(agetodie))
    # print("Day of Death:\n {} \n".format(DOOD))
    # print("Infection array:\n {} \n".format(aliving_infections))
    # print("Total Infections:\n {} \n".format(mic_drop))
    # print("ScareTest:\n {} \n".format(scaretest))

    return aliving_infections, scaretest


# End of defining functions #


# A matrix to record the total value across the ten generations: each row represents a simulation and the four columns are four different output:
# The four output: number of infectious mosquitoes from the wild and released mosquito populations, number of human infections from wild and released mosquitoes
blank_results = np.zeros((N_simulation, 4))  # first number equals size of Z below
stochastic_array = np.asarray(blank_results)

# A matrix to record the per-generation results
Blank_Per_Gen = np.zeros((N_simulation, 40))  # first number equals size of Z below
GenerationMatrix = np.asarray(Blank_Per_Gen)

z = 0

while z < N_simulation:  # Z is the number of simluations

    for J in range(0, 10):  # the 10 generations
        # extract the proportional VC of the wild mosquitoes in that generation
        VC_genX_factor = VC_decline[J]

        if J > 0:  # makes HI proportion of gen 1 infections
            hi = 0.05 * ((GenerationMatrix[z, J - 1] +
                          GenerationMatrix[z, J + 10 - 1]) / GenerationMatrix[z, 0])
        else:
            hi = 0.05

        print(hi)

        x = 0  # wild population loop
        output_array = np.zeros(10)
        dangercount = np.zeros(1)
        while x < N_mosquito:  # adjust this for number of wild mosquitos
            # Where the VC generation factor is applied
            wild = one_mosquito(hi, 0, vcw * VC_genX_factor, sdVCw)
            output_array += wild[0]  # the number of human infections in the 10 gonotrophic cycle
            dangercount += wild[1]  # the number of infectious mosquitoes
            x += 1
        wild_mic_drop = np.sum(output_array)  # total number of human infections by all mosquitoes
        stochastic_array[z, 0] += wild_mic_drop
        stochastic_array[z, 1] += dangercount

        xx = 0  # released population loop
        output_array2 = np.zeros(10)
        releasedanger = np.zeros(1)
        while xx < (N_release * releaseratio):  # adjust this for number of released mosquitos
            if xx < (N_release * releaseratio * fed):
                # the released mosquitoes that are already fed in the lab
                released = one_mosquito(hi, 1, vcr, sdVCr)
            else:
                released = one_mosquito(hi, 0, vcr, sdVCr)
                # the released mosquitoes that faild to fed in the lab, which is roughly equal to wild mosquitoes that have lower vector competence
            output_array2 += released[0]
            releasedanger += released[1]
            xx += 1
        released_mic_drop = np.sum(output_array2)
        stochastic_array[z, 2] += released_mic_drop
        stochastic_array[z, 3] += releasedanger

        GenerationMatrix[z, J] = wild_mic_drop
        GenerationMatrix[
            z, J + 10] = released_mic_drop  # Sets up a 40 column output, first 10 for wild , last 10 released
        GenerationMatrix[z, J + 20] = dangercount
        GenerationMatrix[z, J + 30] = releasedanger

    # print("Total Infection Array, WILD: \n {} \n".format(output_array))
    print("Total Infections WILD:\n {}".format(wild_mic_drop))
    # print("Total Infection Array, RELEASED: \n {} \n".format(output_array2))
    print("Total Infections RELEASED:\n {}".format(released_mic_drop))

    print(GenerationMatrix[z, ])
    print('\n')

    z += 1

print("Infections per generation:\n {}".format(GenerationMatrix))

# print("Wild, Released:\n {} {}".format(stochastic_array[0], stochastic_array[2]))

np.savetxt(f'simulation_BloodfedFemales_Total_{name:03.0f}.csv', stochastic_array, delimiter=",")
np.savetxt(f'simulation_BloodfedFemales_Gens_{name:03.0f}.csv', GenerationMatrix, delimiter=",")

# pads file names with 3 zeros so that it naturally sorts in order (rather than 1,10,2,20 etc)


# hist=plt.hist(stochastic_array[:,[0,2]], 150, histtype = 'bar')
# labels= ["Wild", "Released"]
# plt.ylabel("Frequency")
# plt.xlabel("Infections")
# plt.legend(labels= ["wild","Released"])
# plt.savefig(f'Figure1_{name:03.0f}.png')

# plt.figure()
# hist1=plt.hist(stochastic_array[:,2], 50, histtype = 'bar', label=["Released"], color='orange')
# plt.ylabel("Frequency")
# plt.xlabel("Infections")
# plt.legend(labels= ["Released"])
# plt.savefig(f'Figure2_{name:03.0f}.png')

# plt.figure()
# hist2=plt.hist(stochastic_array[:,[1,3]], 150, histtype = 'bar')
# labels2= ["Wild", "Released"]
# plt.ylabel("Frequency")
# plt.xlabel("Number of Infectious Mosquitoes")
# plt.legend(labels= ["wild","Released"])
# plt.savefig(f'Figure3_{name:03.0f}.png')

print("--- %s seconds ---" % (time.time() - start_time))
