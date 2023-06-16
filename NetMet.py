# -*- coding: utf-8 -*-
import pandas as pd
import pickle
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import datetime, os


def is_consist(LeftSides, RightSides, directions, BucketOfCompounds, Reaction_Direction):
    IsConsist = []
    new = BucketOfCompounds
    for row in range(len(LeftSides)):

        if ((set(LeftSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '>'))) or \
                ((set(LeftSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '='))):

            IsConsist.append(1)
            new = new + RightSides[row]
            # print(str(LeftSides[row])+" IS IN "+str(new))

        elif ((set(RightSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '='))):

            IsConsist.append(1)
            new = new + LeftSides[row]
            # print(str(LeftSides[row])+" IS IN "+str(new))

        else:
            IsConsist.append(0)
            # print(str(LeftSides[row])+" is NOT in "+str(new))

    return IsConsist, new


def plot_complimentary(df, df_steps_, input1, input2, outputfolder, analysisMode_, excludeList=False, DB=None):
    def check_comp(l, r, both):
        if (l == 0) and (r == 0) and (both == 0):
            return 0
        elif (l == 0) and (r == 0) and (both == 1):
            return 1
        elif (l == 1) and (r == 0) and (both == 1):
            return 2
        elif (l == 0) and (r == 1) and (both == 1):
            return 2
        elif (l == 1) and (r == 1) and (both == 1):
            return 3

            # Open log file

    file_object = open(outputfolder + 'Simulation_log.txt', 'a')

    flatui = ['white', 'pink', 'grey', 'brown', 'lime', 'silver']
    group_colors = ['green', 'cyan', 'blue', 'lightblue', 'red', 'orange', 'yellow', 'coral']
    group_colors_dict = {}
    for i, vali in enumerate(flatui + group_colors):
        group_colors_dict[i] = vali

    # envs = list(input2["Env_Content"].values)
    # envs = [i.strip().split(" ")[0:] for i in envs]

    envs = input2[0].values.tolist()

    env_names = [i.split(" ")[0] for i in envs]
    envs = [i.strip().split(" ")[1:] for i in envs]

    ## PATCH - FIX INPUT2##
    D = {}
    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        D[str(i[1]).strip()] = str(i[0]).strip()

    for i, vali in enumerate(envs):
        for j, valj in enumerate(vali):
            envs[i][j] = D[str(valj).strip()]

    #####

    file_object.write('####Plot complimentary #### \n')
    file_object.write('####Envs #### \n')
    file_object.write(str(envs) + ' \n')

    organisms = input1[0].values.tolist()
    organisms = [i.split(" ")[0] for i in organisms]
    enzymes = input1[0].values.tolist()
    enzymes = [i.split(" ")[1:] for i in enzymes]

    file_object.write('####organisms #### \n')
    file_object.write(str(organisms) + ' \n')

    file_object.write('####Mode - interaction #### \n')
    file_object.write(str(envs) + ' \n')

    file_object.write('####enzymes #### \n')
    file_object.write(str(enzymes) + ' \n')

    # Dictionary - Compound index to Compound Code
    CompoundDict = {}

    for i in DB['compound_labels_jun.txt'][["Index", "Compound"]].values.tolist():
        CompoundDict[int(i[0])] = i[1]

    for i in DB['compound_labels_jun.txt'][["Compound_Code", "Compound"]].values.tolist():
        CompoundDict[i[0]] = i[1]

    AllCompDict = {}
    for organism in list(df.columns):
        listofdfs = []
        allCols = []
        for env in range(len(df[organism])):
            cols = [CompoundDict[i] for i in df[organism][env]]
            allCols += cols
            allCols = list(set(allCols))

        for env in range(len(df[organism])):
            cols = [CompoundDict[i] for i in df[organism][env]]
            df_comp = pd.DataFrame([[0] * len(allCols)], columns=allCols)
            df_comp[cols] = int(1)
            listofdfs.append(df_comp)
            # make columns of all df_comp identicle before concat

        AllCompDict[organism] = pd.concat(listofdfs).fillna(0).reset_index(drop=True)

    for count_, k in enumerate(envs):
        cmv = [int(vec) for vec in DB['biomass_vector.txt']["biomass_vector"].values.tolist()[0].strip().split(" ")]
        cmv = list(set([CompoundDict[vec] for vec in cmv]))

        listofdfs = []
        for i in AllCompDict.keys():
            tempdf = pd.DataFrame(AllCompDict[i].T[count_]).reset_index()
            tempdf.columns = ["index", i]
            listofdfs.append(tempdf.set_index("index"))

        df_forJoin = pd.DataFrame(index=cmv)
        df_forJoin[0] = 0

        for i, vali in enumerate(listofdfs):
            df_forJoin = df_forJoin.join(vali)

        listofdfs = df_forJoin.fillna(0)
        del listofdfs[0]

        if excludeList != False:
            listofdfs = listofdfs[~listofdfs.index.isin(excludeList)]

        # Check complementarity for each pair.

        if analysisMode_ == "interaction":
            cols = list(listofdfs.columns)
            for i, vali in enumerate(organisms):
                for j, valj in enumerate(organisms[i:]):
                    if vali != valj:
                        listofdfs[vali + "-" + valj] = [check_comp(xx[0], xx[1], xx[2]) for xx in
                                                        listofdfs[[vali, valj, vali + "-" + valj]].values.tolist()]

        # 1 - arrange compounds according to group
        # 2 - add colored column according to group color
        ####
        groups = list(DB["color_index.txt"]["Group"].unique())
        listofdfs = pd.concat([listofdfs, DB["color_index.txt"][["Group", "Compound"]].set_index("Compound")],
                              join="inner", axis=1)
        listofdfs = listofdfs.sort_values(by="Group")
        tempgroup = listofdfs["Group"].values.tolist()
        listofdfs[listofdfs.index.isin([CompoundDict[int(compk)] for compk in k])] = 4
        m = listofdfs[organisms] == 1

        listofdfs[organisms] = listofdfs[organisms].where(~m, 5)
        listofdfs["Group"] = tempgroup

        groupToValDict = {}

        groupToValDict['Not produced'] = 0
        groupToValDict['Complementary'] = 1
        groupToValDict['Produced by one of the species'] = 2
        groupToValDict['Produced by both species'] = 3
        groupToValDict['Environment'] = 4
        groupToValDict['Produced individually'] = 5
        r = list(range(6, 15))
        for i, vali in enumerate(groups):
            groupToValDict[vali] = r[i]

        listofdfs.replace(groupToValDict, inplace=True)

        ####
        fig = plt.figure(figsize=(40, 40))
        ax = plt.subplot(111)

        sns.set(font_scale=3)

        # Change Group with Category
        listofdfs["Category"] = listofdfs["Group"]
        del listofdfs["Group"]

        img = sns.heatmap(listofdfs,
                          square=True,
                          cmap=sns.color_palette(sns.xkcd_palette(flatui + group_colors)),
                          xticklabels=True,
                          yticklabels=True,
                          linewidths=0.1,
                          linecolor='gray',
                          cbar=False)

        img.set_xticklabels(img.get_xticklabels(), rotation=90, fontsize=20)
        img.set_yticklabels(img.get_yticklabels(), fontsize=20)

        img.tick_params(labelsize=20)

        cax = plt.gcf().axes[-1]
        cax.tick_params(labelsize=20)

        # LEGEND
        import matplotlib.patches as mpatches

        colors = flatui + group_colors
        texts = list(groupToValDict.keys())

        patches = [mpatches.Patch(color=colors[:len(flatui)][i], label="{:s}".format(texts[:len(flatui)][i])) for i in
                   range(len(texts[:len(flatui)]))]

        l1 = plt.legend(handles=patches, fancybox=True, facecolor="white", bbox_to_anchor=(1.05, 0.6), loc='lower left',
                        ncol=1, title="Production mode")

        patches = [mpatches.Patch(color=colors[len(flatui):][i], label="{:s}".format(texts[len(flatui):][i])) for i in
                   range(len(texts[len(flatui):]))]
        plt.legend(handles=patches, fancybox=True, facecolor="white", bbox_to_anchor=(1.05, 1.0), loc='upper left',
                   ncol=1, title="Biochemical category of \nessential cellular compounds")

        plt.gca().add_artist(l1)
        fig.savefig(outputfolder + str(env_names[count_]) + ".png", bbox_inches='tight')
        file_object.close()


def simulation(input1, input2, resfolder, analysisMode_, DB):
    # Open log file
    file_object = open(resfolder + 'Simulation_log.txt', 'a')

    organisms = input1[0].values.tolist()
    organisms = [i.split(" ")[0] for i in organisms]
    enzymes = input1[0].values.tolist()
    enzymes = [i.split(" ")[1:] for i in enzymes]

    # Create all pairs##############################
    # Present multiple analysis or a single analysis for each env
    if (analysisMode_ == "interaction"):
        orgs = []
        enzs = []

        for i, vali in enumerate(organisms):
            for j, valj in enumerate(organisms):
                if (i < len(organisms)) and (j > i):
                    # print(vali+" "+valj)
                    orgs.append(vali + "-" + valj)
                    enzs.append(list(set(enzymes[i] + enzymes[j])))

        organisms = organisms + orgs
        enzymes = enzymes + enzs

    file_object.write('#### List of organisms #### \n')
    file_object.write(str(organisms) + ' \n')
    file_object.write('#### List of enzymes #### \n')
    file_object.write(str(enzymes) + ' \n')

    input1_ = pd.DataFrame()
    input1_["Organism"] = organisms
    input1_["Enzymes_List"] = enzymes

    # envs = list(input2["Env_Content"].values)
    # envs = [i.strip().split(" ")[1:] for i in envs]

    envs = input2[0].values.tolist()

    env_names = [i.split(" ")[0] for i in envs]
    envs = [i.strip().split(" ")[1:] for i in envs]

    file_object.write('#### List of Env #### \n')
    file_object.write(str(envs) + ' \n')

    ## PATCH - FIX INPUT2##
    D = {}
    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        D[str(i[1]).strip()] = str(i[0]).strip()

    for i, vali in enumerate(envs):
        for j, valj in enumerate(vali):
            envs[i][j] = D[str(valj).strip()]

    #####

    # Dictionary - relevant enzymes of each organism
    OrgEnzDict = {}

    # Dictionary - Compound index to Compound Code
    CompoundDict = {}

    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        CompoundDict[i[0]] = i[1]

    # enviornments:
    listOfEnvsOfOrganisms = []
    newOrganismsCompounds = []

    SimSteps_df = pd.DataFrame(columns=["Cycle", "Organism", "InitialEnv", "numberOfCompounds", "Compounds"])
    simulation_steps = []

    for count_, k in enumerate(envs):
        file_object.write('#### Simulation - env ' + str(k) + ' \n')
        file_object.write(str(k) + ' \n')

        # For each Organism:
        Enzymes_List = []
        k = [int(kk) for kk in k]
        newOrganismsCompounds = []
        for i, vali in enumerate(input1_["Organism"].values):
            Enzymes_List = input1_[input1_["Organism"] == vali]
            Enzymes_List = list(Enzymes_List["Enzymes_List"].values)[0]
            Enzymes_List = [j for j in Enzymes_List if str(j) != 'nan']

            file_object.write('#### Simulation - Organism ' + str(vali) + ' \n')
            file_object.write(str(vali) + ' \n')

            # Enzymes_list code to index:
            ll = DB["full_enzymes_labels_jun.txt"][["Index", "Enzyme_Code"]].values.tolist()
            d = {}
            for j in ll:
                d[j[1]] = j[0]

            Enzymes_List_ = []
            for j in Enzymes_List:
                try:
                    Enzymes_List_.append(d[j])
                except:
                    # print(str(j)+" Not in list")
                    file_object.write('#### Simulation - Enzyme list : ' + str(j) + ' Not in list \n')

            # Find reactions consist input Compounds: check if all compounds in envirnment (k) are in the same side of the equation
            # AND an enzyme (j) is in the list --> TAKE the other side of equation according to equilibrium > < =

            df = DB["ec_reac_mapping_jun.txt"][DB["ec_reac_mapping_jun.txt"]['Index'].isin(Enzymes_List_)]

            RelevantReactions = df["Reactions"].values.tolist()
            RelevantReactions = [x.lstrip().rstrip().split(":")[1].lstrip().rstrip().split(" ") for x in
                                 RelevantReactions]
            RelevantReactions = [item for sublist in RelevantReactions for item in sublist]

            for x, valx in enumerate(RelevantReactions):
                try:
                    RelevantReactions[x] = int(valx)
                except:
                    # print("cant make "+valx+" int")
                    file_object.write('#### Simulation - Relevant Reactions : ' + 'cant make ' + str(valx) + ' int  \n')

            RelevantReactions = [x for x in RelevantReactions if x != '']

            OrgEnzDict[i] = RelevantReactions

            df = DB["reactions_3_balanced.txt"][
                DB["reactions_3_balanced.txt"]['Reaction_index'].isin(RelevantReactions)]

            RL = df["Reaction_Left"].values.tolist()
            RR = df["Reaction_Right"].values.tolist()
            Reaction_Direction_ = df["Reaction_Direction"].values.tolist()

            OnlyNewCompounds = []
            newCompounds = []
            prevCompounds = k
            C = 0

            IC, newCompounds = is_consist(RL, RR, Reaction_Direction_, prevCompounds,
                                          Reaction_Direction=Reaction_Direction_)
            simulation_steps.append([count_, C, vali, k, len(set(newCompounds)), list(set(newCompounds))])

            OnlyNewCompounds = list(set(newCompounds).intersection(set(prevCompounds)))

            while set(newCompounds) != set(prevCompounds):
                OnlyNewCompounds = list(set(newCompounds).intersection(set(prevCompounds)))
                simulation_steps.append([count_, C, vali, k, len(set(newCompounds)), list(set(newCompounds))])
                prevCompounds = list(set(newCompounds))
                IC, newCompounds = is_consist(RL, RR, Reaction_Direction_, prevCompounds,
                                              Reaction_Direction=Reaction_Direction_)
                C += 1

            file_object.write('#### Simulation : ' + str(C) + ' Cycles \n')

            # Here build final network snapshot (input cytoscape)
            # 1 - Filter by relevant enzyme (done) 2 - filter by relevant reactions (only reactions which are subsets of the final env)
            # 3 - flatten the list

            newCompoundsCodes = [CompoundDict[k] for k in newCompounds]

            if C > 0:
                newOrganismsCompounds.append(list(set(newCompoundsCodes)))
            else:
                newOrganismsCompounds.append([])

        listOfEnvsOfOrganisms.append(newOrganismsCompounds)

    SimSteps_df = pd.DataFrame(simulation_steps,
                               columns=["EnvIndex", "Cycle", "Organism", "InitialEnv", "numberOfCompounds",
                                        "Compounds"])
    Final_df = pd.DataFrame(listOfEnvsOfOrganisms, columns=organisms)

    # Close the log file
    file_object.close()
    return Final_df, SimSteps_df


try:
    BaseFolder = sys.argv[1]
    organismsFileName = sys.argv[2]
    envsFileName = sys.argv[3]
    analysisMode = sys.argv[4]
    UserInputFiles = [organismsFileName, envsFileName]
    # print(sys.argv)
    # print(UserInputFiles)
except:
    print("Error in executing. Using default parameters")
    BaseFolder = "./"
    UserInputFiles = [BaseFolder + "genomes.txt", BaseFolder + "envs.txt"]
    analysisMode = "interaction"
try:
    with open(BaseFolder + "database/DB.pickle", 'rb') as handle:
         DB = pickle.load(handle)
except:
    DB = pd.read_pickle(BaseFolder + "database/DB.pickle")

try:
   os.mkdir(BaseFolder + "Results")
except:
    print()

d = datetime.datetime.today()

FinalFolder = BaseFolder + "Results/" + str(d.year) + str(d.month) + str(d.day) + str(d.hour) + str(d.minute) + str(
    d.second) + "/"
os.mkdir(FinalFolder)

 # 1 - Take Input2 compounds (medium)
Input1 = pd.read_csv(UserInputFiles[0], header=None)
Input2 = pd.read_csv(UserInputFiles[1], header=None)  # , names=["Index","Env_Content"])

# Run simulations
df_final, df_simsteps = simulation(input1=Input1, input2=Input2, resfolder=FinalFolder, analysisMode_=analysisMode)

# Save results to file
df_final.to_csv(FinalFolder + "results.csv")

# df_simsteps.to_csv(FinalFolder+"sim_steps.csv")

# Plot complementarity results and save as png
plot_complimentary(df=df_final,
                   df_steps_=df_simsteps,
                   input1=Input1,
                   input2=Input2,
                   outputfolder=FinalFolder,
                   analysisMode_=analysisMode
                   )

"""
Example execution - 
python3 NetCMPT_v1_15.py ./ ./genomes.txt ./Simple_env.txt single-species
python3 NetCMPT_v1_15.py ./ ./genomes.txt ./Simple_env.txt interaction

"""
