import mitty.lib.generations as gen


def make_generations(gens=1):
  """Give us a mating pair that have common parents only 'gens' generations back."""
  g = [[gen.Sample(0, 0), gen.Sample(0, 1)]] + [[gen.Sample(i + 1, n) for n in range(2 ** (gens - m))] for i, m in enumerate(range(gens))]
  for i, x in enumerate(g[1]):
    x.parents = [g[0][0], g[0][1]]
  for m in range(2, gens + 1):
    for i, x in enumerate(g[m]):
      x.parents = [g[m - 1][2 * i], g[m - 1][2 * i + 1]]
  return g


def print_ancestors(s):
  m = [s]
  while len(m):
    m = [ll for l in m for ll in l.parents]
    print([p.name for p in m])


def test_incest():
  """Stock incest algorithm"""
  # With common parents
  p = [gen.Sample(1, n) for n in range(2)]
  ch = [gen.Sample(2, n) for n in range(2)]
  for n in range(2):
    ch[n].parents = p
  assert gen.incestuous_mating(ch[0], ch[1], incest_generations=1) == True

  # With common grandparents
  gp = [gen.Sample(1, n) for n in range(2)]
  p = [gen.Sample(2, n) for n in range(4)]
  ch = [gen.Sample(3, n) for n in range(8)]

  for n in range(4):
    p[n].parents = gp
  for n in range(8):
    ch[n].parents = p[2*(n/4):2*(n/4) + 2]

  assert not gen.incestuous_mating(ch[0], ch[4], incest_generations=1)
  assert gen.incestuous_mating(ch[0], ch[4], incest_generations=2)


  # With common great-grandparents
  g = make_generations(gens=4)
  assert not gen.incestuous_mating(g[3][0], g[3][1], incest_generations=1)
  assert not gen.incestuous_mating(g[3][0], g[3][1], incest_generations=2)
  assert gen.incestuous_mating(g[3][0], g[3][1], incest_generations=3)


def test_mating():
  """Stock mate choice algorithm."""
  # Create a parent list, then a child list such that there are four pairs, alternately unrelated and related
  p = [gen.Sample(0, n) for n in range(12)]
  c = [gen.Sample(1, n) for n in range(8)]
  c[0].parents = [p[0], p[1]]
  c[1].parents = [p[2], p[3]]  # Unrelated pair
  c[2].parents = [p[4], p[5]]
  c[3].parents = [p[4], p[5]]  # Related pair
  c[4].parents = [p[6], p[7]]
  c[5].parents = [p[8], p[9]]  # Unrelated pair
  c[6].parents = [p[10], p[11]]
  c[7].parents = [p[10], p[11]]  # Related pair

  mating_list, incest_list = gen.mate_in_sequence(list(c), incest_generations=1)
  assert mating_list == [[c[5], c[4]], [c[1], c[0]]], mating_list
  assert incest_list == [[c[7], c[3]], [c[6], c[2]]], incest_list
