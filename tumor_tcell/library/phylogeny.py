import itertools

def get_phylogeny(agent_ids):
  # agent_ids: list of string values
  # make phylogeny with {mother_id: [daughter_1_id, daughter_2_id]}
  phylogeny = {agent_id: [] for agent_id in agent_ids}
  for agent1, agent2 in itertools.combinations(agent_ids, 2):
      if agent1 == agent2[0:-1]:
          phylogeny[agent1].append(agent2)
      elif agent2 == agent1[0:-1]:
          phylogeny[agent2].append(agent1)
  return phylogeny


def daughter_ab(mother_id):
    return [
        str(mother_id) + "A",
        str(mother_id) + "B"
    ]
