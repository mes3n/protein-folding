import matplotlib.pyplot as plt
import json

with open('results/raw.json', 'r') as f:
  data = json.load(f)

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

ref_similarity = data['normal i=500']['ref_similarity']
step_similarity = data['normal i=500']['step_similarity']

plt.plot(step_similarity)
plt.ylim(0, max(step_similarity) * 1.2)
plt.xlabel('Iterations')
plt.ylabel('Distance [Å]')
plt.title('Average distance of each atom to previous iteration')
plt.savefig('results/step.svg')
plt.show()
plt.close()

plt.plot(ref_similarity)
plt.ylim(0, max(ref_similarity) * 1.2)
plt.xlabel('Iterations')
plt.ylabel('Distance [Å]')
plt.title('Average distance of each atom to reference protein')
plt.savefig('results/reference.svg') 
plt.show()
plt.close()
