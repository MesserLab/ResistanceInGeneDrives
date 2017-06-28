#!/usr/bin/env python

def parse(out):
    lines = out.split("\n")
    driver = []
    wild = []
    r01 = []
    r10 = []
    r11 = []
    for line in lines:
        if line.startswith("#OUTPUT:"):
            freqs = line.split()
            driver.append(freqs[1])
            wild.append(freqs[2])
            r01.append(freqs[3])
            r10.append(freqs[4])
            r11.append(freqs[5])
    if len(driver) < 40:
        driver + driver[-1]*(40-len(driver))
        wild + wild[-1]*(40-len(driver))
        r01 + r01[-1]*(40-len(driver))
        r10 + r10[-1]*(40-len(driver))
        r11 + r11[-1]*(40-len(driver))
    return [driver, wild, r01, r10, r11]

def main():
    print("test")

if __name__=="__main__":
    main()
