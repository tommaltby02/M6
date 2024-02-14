import numpy as np

analytical = 1

def TrapezoidalIntegration(no_samples, a, b):
        h = (b - a) / no_samples

        result = 0.5 * (F(a) + F(b))
        for i in range(1, no_samples):
                result += F(a + i*h)
        result*= h
        return result
        

def MCIntegration(weighting_function, no_samples):
        samples = np.random.uniform(0, 1, no_samples)
        values = F(samples)
        weights = GetWeightingFunction(weighting_function, samples)

        result = np.mean(values / weights)
        return result

def GetWeightingFunction(type, samples):
        if type == "1":
                return 1
        elif type == "2":
                return 2 * samples
        elif type == "3":
                return 4 * samples**3 
        elif type == "4":
               return 3 * samples**2
def F(x):
        return 3 * x**2

def Run(type, weighting_function, no_samples, no_runs):
        results = []
        for i in range(no_runs):
            if type == "MC":
                results.append(MCIntegration(weighting_function, no_samples)) 
            elif type == "trapezoid":
                results.append(TrapezoidalIntegration(no_samples, 0, 1))
        mean = np.mean(results)
        std = np.std(results)
        mae = np.mean(np.abs(results - np.array(analytical)))
        print(f"Mean integral = {mean}, standard deviation = {std}, MAE = {mae}")
        return

Run("MC", "4", 1000, 100)