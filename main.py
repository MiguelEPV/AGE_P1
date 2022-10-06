import requests
import random
import matplotlib
import time
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
website = "http://memento.evannai.inf.uc3m.es/age/alfa?c="


def dec_to_bin(num):
    return bin(num).replace("0b", "")


def lowest_chrom():
    lowest = 999999999
    chrom_num = 0
    for i in range(pow(2, 32)):
        chromosome = dec_to_bin(i)
        contents = requests.get(website + chromosome)
        res = float(contents.text)
        if res < lowest:
            lowest = res
            chrom_num = i
    print("Valor obtenido " + str(res))
    print("Chromosome: " + str(chrom_num))


def create_individual(size, lst):
    if size == 0:
        return lst
    else:
        lst.append(random.randint(0, 1))
        return create_individual(size - 1, lst)


def init_vectores(pob, n):
    vec = [create_individual(n, []) for x in range(0,pob)]
    return vec

def evaluar_pob(pob):
    values = []
    for ind in pob:
        chrom = ''.join(map(str, ind))
        contents = requests.get(website + chrom)
        res = float(contents.text)
        values.append(res)
    return values


def torneo(pob, participants):
    winners = []
    while len(winners) < len(pob):
        indvs = []
        id_lst = []
        for x in range(0,participants):
            id = random.randint(0,(len(pob)-1))
            id_lst.append(id)
            indvs.append(pob[id])
        win = min(indvs)
        win_id = indvs.index(win)
        winners.append(id_lst[win_id])

    return winners


def cruzar(parents, genes, cross_prob, mut_prob):
    new_pob = []
    for i in range(0, len(parents), 2):
        dad_id = parents[i]
        mom_id = parents[i+1]
        dad = genes[dad_id]
        mom = genes[mom_id]
        child_x = []
        child_y = []
        for j in range(0, len(dad)):
            rand_x = random.random()
            rand_y = random.random()
            mut_x = random.random()
            mut_y = random.random()
            if rand_x <= cross_prob:
                child_x.append(dad[j])
            else:
                child_x.append(mom[j])

            if rand_y <= cross_prob:
                child_y.append(mom[j])
            else:
                child_y.append(dad[j])

            if mut_x <= mut_prob:
                if child_x[j] == 1:
                    child_x[j] = 0
                else:
                    child_x[j] = 1

            if mut_y <= mut_prob:
                if child_y[j] == 1:
                    child_y[j] = 0
                else:
                    child_y[j] = 1
        new_pob.append(child_x)
        new_pob.append(child_y)
    return new_pob


def cruzar_2(parents, genes, partir, cross_prob, mut_prob):
    new_pob = []
    for i in range(0, len(parents), 2):
        dad_id = parents[i]
        mom_id = parents[i + 1]
        dad = genes[dad_id]
        mom = genes[mom_id]
        child_x = []
        child_y = []
        # cruce = 0
        count = 1
        sections = len(dad)/partir
        rand_num1 = random.random()
        rand_num2 = random.random()
        for j in range(0, len(dad)):
            mut_x = random.random()
            mut_y = random.random()
            if rand_num1 <= cross_prob:
                child_x.append(dad[j])
            else:
                child_x.append(mom[j])

            if rand_num2 <= cross_prob:
                child_y.append(mom[j])
            else:
                child_y.append(dad[j])
            if count % sections == 0:
                count = 0
                rand_num1 = random.random()
                rand_num2 = random.random()
                slice_start = j - (int(sections)-1)
                slice_end = j+1
                if mut_x <= mut_prob:
                    child_x[slice_start:slice_end] = mutar_2(child_x[slice_start:slice_end])
                if mut_y <= mut_prob:
                    child_y[slice_start:slice_end] = mutar_2(child_y[slice_start:slice_end])
            # elif count%sections == 0 and cruce == 1:
            #     cruce = 0
            #     count = 0

            count +=1

            # if mut_x <= mut_prob:
            #     if child_x[j] == 1:
            #         child_x[j] = 0
            #     else:
            #         child_x[j] = 1
            #
            # elif mut_y <= mut_prob:
            #     if child_y[j] == 1:
            #         child_y[j] = 0
            #     else:
            #         child_y[j] = 1
        new_pob.append(child_x)
        new_pob.append(child_y)
    return new_pob

def mutar_2(segment):
    chrom = ''.join(map(str, segment))
    chrom = '0b' + chrom
    dec_value = int(chrom, 2)
    # mut_num = int(dec_value/2)
    # # mut_num = int(dec_value*0.1)
    mut_num = 25
    rand_x = random.randint(0,mut_num)
    new_segment = []
    prob = random.randint(0,1)

    if prob == 0:
        dec_value = (dec_value - rand_x) % (2**len(segment))
    else:
        dec_value = (dec_value + rand_x) % (2**len(segment))
    new_genes = format(dec_value, "b").zfill(len(segment))
    for i in new_genes:
        new_segment.append(int(i))
    return new_segment



def main():
    pob_num = 200
    genes = 80
    generations = 50
    segments = 10
    mutation_prob = 0.2
    cross_prob = 0.5
    percentage = 3/100
    participants = int(pob_num * percentage)
    graph_x = []
    graph_y = []
    pob = init_vectores(pob_num, genes)
    # with open('data.txt', 'w') as f:
    #     f.write("Poblacion: " + str(pob_num) + ", Generaciones: " + str(generations) +
    #             ", Probabilidad de Cruce: " + str(cross_prob) + ", Probabilidad de Mutacion: " + str(mutation_prob))
    #     f.write('\n')
    #     f.write(str(pob))
    #     f.write('\n')
    #     f.write('----------------------------------------------')
    #     f.write('\n')
    for i in range (0, generations):
        t1 = time.time()
        graph_x.append(i)
        eval = evaluar_pob(pob)
        graph_y.append(min(eval))
        best = torneo(eval, participants)
        cross = cruzar_2(best, pob, segments, cross_prob, mutation_prob)
        # cross = cruzar(best, pob, cross_prob, mutation_prob)
        pob = cross
        print("Gen: " + str(i) + " Best value: " + str(min(eval)) + " Mean: " + str(np.mean(eval)))
        # if round(np.mean(eval),2) == round(min(eval), 2):
        #     mutation_prob += 0.02
        t2 = time.time() - t1
        print("Time: " + str(t2))
    final_eval = evaluar_pob(pob)
    print(min(final_eval))
    plt.plot(graph_x, graph_y)
    plt.xlabel('Generation')
    plt.ylabel('Best Value')
    plt.show()


if __name__ == "__main__":
    main()










