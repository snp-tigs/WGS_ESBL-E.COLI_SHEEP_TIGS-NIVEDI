import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

#\data as float values
data_series = [
    #rectal_swab
    [
2.950019884,
2.90276629,
2.907175548,
3.011468462,
3.012839934,
2.860918694,
2.897405604,
2.944364264,
3.060514332,
2.951025973,
2.952001091,
2.978441714,
2.9601435,
3.039266524,
2.920460765,
3.01614998,
2.875146633,
3.016472124,
2.935232469,
2.913452632,
2.876865664,
2.973606292,
3.061878784,
2.898610751,
2.987110408,
3.12415766,
3.0010172,
2.951765536,
3.023514998,
3.021309279,
2.968427002,
3.083411036,
2.927796385,
2.94911442,
2.964586303,
2.883278455,
2.984009766,
2.964324161,
2.945332107,
2.948016532,
3.020767683,
2.924989407,
3.02513302,
2.944744897,
2.892026979,
2.875006052,
2.971038907,
2.746747206,
3.005758008,
2.763066302,
2.996852039,
2.508593134,
2.335026393,
2.128079507,
1.353975108],
#env
    [
2.924139506,
2.925151769,
3.017537389,
2.956829703,
2.940178389,
3.018314998,
2.946547922,
2.958594478,
2.967616978,
2.43671079,
2.864116428,
2.989267064,
2.976481248,
2.932876604,
2.858105017,
2.89719468,
1.769803198]
]

series_names = ['Sheep_isolates', 'Environmental_isolates']
box_colors = ['#EDCBA8', '#B8D1AE']

plt.figure(figsize=(8, 6))
sns.boxplot(data=data_series, palette=box_colors)
sns.stripplot(data=data_series, color='black', jitter=True, size=5)

plt.xlabel('sample type')
plt.ylabel('Shannon diversity drug classes based on RPM')
plt.xticks(ticks=range(len(series_names)), labels=series_names)

def add_significance_annotation(ax, x1, x2, y, significance):
    h = 0.003
    col = 'k'
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
    ax.text((x1 + x2) * 0.5, y + h, significance, ha='center', va='bottom', color=col)

ax = plt.gca()
y_max = np.max([np.max(series) for series in data_series])
comparisons = [(0, 1)]

# Perform Kruskal-Wallis H-test
# *data_series unpacks the list of lists into individual arguments
h_stat, p_val = stats.kruskal(*data_series)

# Determine significance level
if p_val <= 0.001:
    significance = '***'
elif p_val <= 0.01:
    significance = '**'
elif p_val <= 0.05:
    significance = '*'
else:
    significance = 'ns'

# Add significance level annotation
for i, (x1, x2) in enumerate(comparisons):
    y = y_max + (i * 0.01)
    add_significance_annotation(ax, x1, x2, y, significance)

# Output results to screen
print(f"Kruskal-Wallis H-statistic: {h_stat:.4f}")
print(f"P-value: {p_val:.4f}")
print(f"Significance: {significance}")

plt.show()
