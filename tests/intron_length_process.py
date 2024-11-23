from statistics import mean, median
import numpy

if __name__ == "__main__":

    f = open(r"C:\Users\tonyx\Downloads\intron_lengths.txt")
    line = f.readline()

    arr = []

    while line != "":

        if line.strip().isnumeric():
            arr.append(int(line.strip()))

        line = f.readline()

    print("mean: " + str(mean(arr)))
    print("median: " + str(median(arr)))
    print("max: " + str(max(arr)))
    print("min: " + str(min(arr)))

    x = numpy.quantile(arr, [0, 0.25, 0.5, 0.75, 0.95, 1])
    print(x)

    f2 = open(r"C:\Users\tonyx\Downloads\new_introns.txt", "w")
    arr.sort()
    for item in arr:
        f2.write(str(item) + "\n")

    f2.close()


