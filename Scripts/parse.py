#!/usr/bin/env python

def parse(out):
    lines = out.split("\n")
    driver = []
    wild = []
    r1 = []
    r2 = []
    for line in lines:
        if line.startswith("#OUTPUT:"):
            freqs = line.split()
            driver.append(freqs[1])
            wild.append(freqs[-1])
            r1.append(freqs[2])
            if len(line.split()) > 2:
                r2.append(freqs[3])
    if len(driver) < 40:
        driver + driver[-1]*(40-len(driver))
        wild + wild[-1]*(40-len(driver))
        r1 + r1[-1]*(40-len(driver))
        r2 + r2[-1]*(40-len(driver))
    return [driver, wild, r1, r2]

def main():
    print("test")

if __name__=="__main__":
    main()
