import numpy as np
import numexpr as ne
# import math as m
# import sys as sys
# import matplotlib.pyplot as plt
import time
start_time = time.time()

np.set_printoptions(threshold=np.inf)


"""
Section 1: specify simulation scale and import model parameters
"""

N_mosquito = 5000  # Number of wild mosquito populations
N_release = N_mosquito * 0.1  # Number of released mosquitoes
N_simulation = 100  # Number of simulations in each run

modelnumber = 0  # use this to run on laptop
# modelnumber=int(sys.argv[1])          # for cluster

# An matrix with all parameter values: each row is a parameter set and each column is one parameter
parameterarray = np.loadtxt("SensitivityAnalysisBloodfed", delimiter=",")
currentparameters = parameterarray[modelnumber, :]

with np.printoptions(precision=3, suppress=True):
    print(currentparameters)  # Print the current paramter values

hi, FF, vcw, vcr, ovi, bitephase, EIP, mubites, releaseratio, a, b, s, c, sdovi, sdhost, sdEIP, sdVCw, sdVCr, fed, name, VC1, VC2, VC3, VC4, VC5, VC6, VC7, VC8, VC9, VC10 = currentparameters
# Assign parameter values to variables for easy use in the rest of the script
VC_decline = np.array([VC1, VC2, VC3, VC4, VC5, VC6, VC7, VC8, VC9, VC10])

"""
parameter_dict = dict(hi=hi,
                      FF=FF,
                      vcw=vcw,
                      vcr=vcr,
                      ovi=ovi,
                      bitephase=bitephase,
                      EIP=EIP,
                      mubites=mubites,
                      releaseratio=releaseratio,
                      a=a, b=b, s=s, c=c,
                      sdovi=sdovi,
                      sdhost=sdhost,
                      sdEIP=sdEIP,
                      sdVCw=sdVCw,
                      sdVCr=sdVCr,
                      fed=fed)
"""

"""
Parameter explanation:
# hi                    = % of Humans Infectious (dafault = 0.05)
# FF                  = age of the mosquitoes when released (default = 5, i.e. the released mosquitoes are five days old)
# vcw                = mean vector competence of the wild mosquitoes (default = 0.524)
# vcr                  = mean vector competence of the released mosquitoes (default = 0.2)
# ovi                  = mean duration of the oviposition period (default = 4 days)
# bitephare        = mean duration of the host-seeking period (default = 4 days)
# EIP                 = mean duration of the extrinsic incubation period (the time between a mosquito acquires the virus and it becomes infectious, default = 11 days)
# mubites           = mean number of bites in each host-seeking period (default = 2)
# releaseratio     = size of the released population relative to the wild population (default = 1, which represents releasing 10% of the wild population size)
# a, b, s, c           = parameters that define the mosquito age-dependent survival function
# sdovi               = standard deviation of the oviposition period (default = 1 days)
# sdhost              = standard deviation of the host-seeking period (default = 1 days)
# sdbite              = standard deviation of the number of bites in each host-seeking period (default = 1)
# sdEIP               = standard deviation of EIP (default = 1 days)
# sdVCw             = standard deviation of the wild mosquito VC (default = 0.125)
# sdVCr              = standard deviation of the released mosquito VC (default = 0.125)
# fed                   = proportion of the released mosquitoes being fed before the release (default = 1)
# name               = index of runs (for tracking different runs on the cluster, not part of the simulation but is used for writing simulation output)
# VC1 - VC10   = the VC of the wild mosquito population throughout the 10 generations, expressed as proportion to the initial wild population VC (i.e. the first element is always 1)

# Note 1: vector competence represents the mosquito's ability to transmit diseases. We assume each mosquitoes have a VC range between 0 and 1.
# Note 2: #    we are proposing to release female mosquitoes after blood-feeding them in the lab, but practically some individuals might fail to feed in the lab

# The change of VC across the 10 generations for four difference scenarios at baseline
# VC_decline=np.array ([1,0.959,0.921,0.886,0.853,0.822,0.794,0.767,0.743,0.719])       # Only males
# For only males, there're effectively 0 released mosq (only VC change)
# VC_decline=np.array ([1,0.928,0.865,0.810,0.761,0.718,0.680,0.646,0.615,0.588])       # Blood-fed Fem + Males
# VC_decline=np.array ([1,0.956,0.916,0.879,0.845,0.813,0.783,0.756,0.730,0.707])       #  Unfed Fem + Males
# VC_decline=np.array ([1,1,1,1,1,1,1,1,1,1])                                                                      # no releases
"""


"""
Section 2: Simulation functions
"""


# Function to calulate the mortality of mosquitoes depending on their age
def age_dependence(x, a1, b1, s1, c1, delta_t):
    """
    x: age of the mosquito
    a1, b1, s1, c1: parameters to define the hazard function (Styler et al. 2007)
        a1 = 1.799999999999999951e-03
        b1 = 1.416000000000000036e-01
        s1 = 1.072999999999999954e+00
        c1 = 2.199999999999999872e-02
        adjust c to change agerage life expectancy, c=0.022 is 21 days
    delta_t: time interval
    """
    return ne.evaluate('exp(-((a1 * exp(b1 * x)) / (1 + (a1 * s1 / b1) * (exp(b1 * x) - 1)) + c1) * delta_t)')


# probability of surviving any 0.1 day interval:
age_wild = age_dependence(np.arange(0, 100, 0.1), a, b, s, c, 0.1)
age_wild[-1] = 0  # force the mosquito to die at the last time breakpoint

# probability of surviving to any 0.1 day for wild mosquitoes:
survival_wild = np.cumprod(age_wild)

age_release = age_wild[int(10*FF):]  # released mosquitoes have no mortality before being released
# probability of surviving to any 0.1 day for released mosquitoes (since the time of release):
survival_release = np.cumprod(age_release)


# Function to sample the life span of a mosquito
def age_of_death(survival_p):  # age of death calculator
    # survival_p: option = survival_wild or survival_release
    rando = np.random.random_sample(1)
    itemindex = np.asarray(np.where(survival_p < rando))
    DEATH = itemindex[0, 0] / 10
    return DEATH


# Function to simulate one single mosquito
# We first simulate a sequence of 10 gonotrophic cycle and them use the time of death to clip the sequence to preserve only the parts when the mosquito is alive
def one_mosquito(hi, trinary, vcmu, vcsigma):
    # trinary: 0=wild, 1=release after successful blood-feeding, 2=release after failed bloodfeeding

    # Vector Competence:
    mu, sigma = vcmu, vcsigma
    vc = np.random.normal(mu, sigma, 1)  # each individual mosquito has a constant VC

    # Duration of Oviposition:
    mu, sigma = ovi, sdovi  # Duration of Oviposition
    # 10 numbers to represent 10 gonotrophic cycles (each includes a host-seeking and an oviposition period)
    do = np.random.normal(mu, sigma, 10)
    # assume that each oviposition period must be at least 3 days and no more than 30 days
    trunc_do = np.clip(do, 3, None)
    ado = np.asarray([trunc_do])

    # For each gonotrophic cycle, we start with the oviposition phase and then the host-seeking phase
    # For wild mosquitoes, it requires a two-day maturation before the first host-seeking, which equilvalent to a two-day oviposition period
    # For consistency with the release mosquitoes, which enters the simulation in the oviposition period, we consider this "two-day" period in wild mosquitoes as its first oviposition peirod
    # So the next few lines change the first oviposition peirod of the wild mosquitoes
    firstO = np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1])  # Activated only for the wild pop
    add_2 = np.array([2, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    ado0 = ado * firstO
    ado2 = np.asarray(ado0 + add_2)

    if trinary == 1:
        final_ato = ado  # release females successfully fed
    elif trinary == 2:
        final_ato = ado0  # release females unsuccessfully fed
    else:
        final_ato = ado2  # wild females

    atdo = final_ato.transpose()
    attdo = atdo.transpose()

    # Duration of Bloodfeeding:
    mu, sigma = bitephase, sdhost
    db = np.random.normal(mu, sigma, 10)
    # At least spending one day in blood-feeding
    trunc_db = np.clip(db, 1, None)
    adb = np.asarray([trunc_db])
    atdb = np.transpose(adb)

    # Bites per gonotrophic cycle (i.e. each host-seeking period):
    bites = np.random.poisson(mubites, 10)
    nb = np.array([bites])
    anb = np.round(nb)  # rounded bc: bites must be integers
    trunc_anb = np.clip(anb, 1, None)  # at least one bites per cycle and at most 10 bites
    atnb = trunc_anb.transpose()

    # Extrinsic Incubation Period:
    mu, sigma = EIP, sdEIP
    eip_draw = np.random.normal(mu, sigma, 1)  # Only need one value for each individual
    eip = np.clip(eip_draw, 5, None)  # EIP at least 5 days and at most 30 days

    # a vector to record whether the mosquito got infected during each gonotrophic cycle
    infection_status = np.zeros([10])

    # cycle through the first nine gonotrophic cycle:
    for x in range(0, 9):
        bitesincyc = int(atnb[x])  # extract the number of bites in this gonotrophic cycle
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

    # Create a counter to record the days since the mosquitoes got affected for each gonotrophic cycle:
    # if counter > EIP, the mosquito becomes infectious
    empty_list = np.zeros(10)
    empty_array = np.array([empty_list])
    counter = empty_array.transpose()

    for x in range(0, 9):  # Calculates initial counter input
        biteday = (tIS[x])
        # calls the infectious binary signal
        bitephaseday = (atdb[x])
        # calls the length of biting cycle when mosquito was infected
        randn = np.random.random_sample(())
        if biteday == 1:
            counter[x] = randn * bitephaseday
            # Within the host-seeking perid in which the mosquito got infected, what is the time between the infectious bites and the end of the biting phase

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

    # Number of humans infected by this mosquito in each gonotrophic cycle:
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
    if trinary == 0:  # Links the dif DOD functions to trinary, trinary = 1 indicates that we are simulating a release mosquito
        DOOD = age_of_death(survival_wild)
    else:
        DOOD = age_of_death(survival_release)
    # selects the death age function, assumes 5 days old for the released

    for x in range(0, 10):  # eliminates all transmissions after death
        biteage = abiting_age[x]  # age of the mosquito at the end of a gonotrophic cycle
        if biteage <= DOOD:  # If the mosquito has not died
            agetodie[x] = 1  # makes a binary check 1=alive
            if x < 9:
                # ensures that x+1 isn't out of range
                # the mosquito enters the next gonotrophic cycle, which means it has the chance to bite more people
                agetodie[x + 1] = 1
    # After this loop, agetodie should indicate whether the mosquito is still alive at the beginning of next gonotrophic cycle,
    # which gives the mosquito chance to bite more humans in the next gonotropic cycle

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


"""
Section 3: simulate mosquito releases and disease transmission
"""


# A matrix to record the total value across the ten generations:
# each row represents a simulation and the four columns are four different output:
# The four output: number of infectious mosquitoes from the wild and released mosquito populations, number of human infections from wild and released mosquitoes
blank_results = np.zeros((N_simulation, 4))  # first number equals size of Z below
stochastic_array = np.asarray(blank_results)

# A matrix to record the per-generation results
Blank_Per_Gen = np.zeros((N_simulation, 40))  # first number equals size of Z below
GenerationMatrix = np.asarray(Blank_Per_Gen)

z = 0

while z < N_simulation:  # Z is the number of simluations
    print("Simulation number", z)

    for J in range(0, 10):  # the 10 generations
        print("Generation", J)

        # extract the proportional VC of the wild mosquitoes in that generation
        VC_genX_factor = VC_decline[J]

        # update human infection proportion: hi
        if J > 0:  # makes HI proportion of gen 1 infections
            hi = 0.05 * ((GenerationMatrix[z, J - 1] +
                          GenerationMatrix[z, J + 10 - 1]) / GenerationMatrix[z, 0])
        else:
            hi = 0.05

        print("human infection proportion:", hi)

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
                # the released mosquitoes that faild to fed in the lab
                released = one_mosquito(hi, 2, vcr, sdVCr)
            output_array2 += released[0]
            releasedanger += released[1]
            xx += 1
        released_mic_drop = np.sum(output_array2)
        stochastic_array[z, 2] += released_mic_drop
        stochastic_array[z, 3] += releasedanger

        GenerationMatrix[z, J] = wild_mic_drop
        # Sets up a 40 column output, first 10 for wild , last 10 released
        GenerationMatrix[z, J + 10] = released_mic_drop
        GenerationMatrix[z, J + 20] = dangercount
        GenerationMatrix[z, J + 30] = releasedanger

    # print("Total Infection Array, WILD: \n {} \n".format(output_array))
    print("Total Infections WILD:\n {}".format(wild_mic_drop))
    # print("Total Infection Array, RELEASED: \n {} \n".format(output_array2))
    print("Total Infections RELEASED:\n {}".format(released_mic_drop))

    print(GenerationMatrix[z, ])
    print('\n\n\n\n\n')

    z += 1

print("Infections per generation:\n {}".format(GenerationMatrix))

# print("Wild, Released:\n {} {}".format(stochastic_array[0], stochastic_array[2]))

np.savetxt(f'simulation_BloodfedFemales_Total_{name:03.0f}.csv', stochastic_array, delimiter=",")
np.savetxt(f'simulation_BloodfedFemales_Gens_{name:03.0f}.csv', GenerationMatrix, delimiter=",")

# pads file names with 3 zeros so that it naturally sorts in order (rather than 1,10,2,20 etc)

"""
hist=plt.hist(stochastic_array[:,[0,2]], 150, histtype = 'bar')
labels= ["Wild", "Released"]
plt.ylabel("Frequency")
plt.xlabel("Infections")
plt.legend(labels= ["wild","Released"])
plt.savefig(f'Figure1_{name:03.0f}.png')

plt.figure()
hist1=plt.hist(stochastic_array[:,2], 50, histtype = 'bar', label=["Released"], color='orange')
plt.ylabel("Frequency")
plt.xlabel("Infections")
plt.legend(labels= ["Released"])
plt.savefig(f'Figure2_{name:03.0f}.png')

plt.figure()
hist2=plt.hist(stochastic_array[:,[1,3]], 150, histtype = 'bar')
labels2= ["Wild", "Released"]
plt.ylabel("Frequency")
plt.xlabel("Number of Infectious Mosquitoes")
plt.legend(labels= ["wild","Released"])
plt.savefig(f'Figure3_{name:03.0f}.png')
"""

print("--- %s seconds ---" % (time.time() - start_time))
