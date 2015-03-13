import mitty.lib.generations as gen


def make_generations(gens=1):
  """Give us a mating pair that have common parents only 'gens' generations back."""
  g = [[gen.Sample(0, 0), gen.Sample(0, 1)]] + [[gen.Sample(i + 1, n) for n in range(2 ** (gens - m))] for i, m in enumerate(range(gens))]
  for i, x in enumerate(g[1]):
    x.p1, x.p2 = g[0][0], g[0][1]
  for m in range(2, gens + 1):
    for i, x in enumerate(g[m]):
      x.p1, x.p2 = g[m - 1][2 * i], g[m - 1][2 * i + 1]
  return g


def print_ancestors(s):
  m = [s]
  while len(m):
    m = [ll for l in m for ll in [l.p1, l.p2]]
    print([p.name for p in m])


def test_incest():
  """Stock incest algorithm"""
  # With common parents
  p = [gen.Sample(1, n) for n in range(2)]
  ch = [gen.Sample(2, n) for n in range(2)]
  for n in range(2):
    ch[n].p1, ch[n].p2 = p[0], p[1]

  mtr = gen.StockMater(incest_generations=1)
  assert mtr.incestuous_mating(ch[0], ch[1]) == True

  # With common grandparents
  gp = [gen.Sample(1, n) for n in range(2)]
  p = [gen.Sample(2, n) for n in range(4)]
  ch = [gen.Sample(3, n) for n in range(8)]

  for n in range(4):
    p[n].p1, p[n].p2 = gp[0], gp[1]
  for n in range(8):
    ch[n].p1, ch[n].p2 = p[2*(n/4)], p[2*(n/4) + 1]

  mtr = gen.StockMater(incest_generations=1)
  assert not mtr.incestuous_mating(ch[0], ch[4])

  mtr = gen.StockMater(incest_generations=2)
  assert mtr.incestuous_mating(ch[0], ch[4])

  # With common great-grandparents
  g = make_generations(gens=4)

  mtr = gen.StockMater(incest_generations=1)
  assert not mtr.incestuous_mating(g[3][0], g[3][1])

  mtr = gen.StockMater(incest_generations=2)
  assert not mtr.incestuous_mating(g[3][0], g[3][1])

  mtr = gen.StockMater(incest_generations=3)
  assert mtr.incestuous_mating(g[3][0], g[3][1])


def test_mating():
  """Stock mate choice algorithm."""
  # Create a parent list, then a child list such that there are four pairs, alternately unrelated and related
  p = [gen.Sample(0, n) for n in range(12)]
  c = [gen.Sample(1, n) for n in range(8)]
  c[0].p1, c[0].p2 = p[0], p[1]
  c[1].p1, c[1].p2 = p[2], p[3]  # Unrelated pair
  c[2].p1, c[2].p2 = p[4], p[5]
  c[3].p1, c[3].p2 = p[4], p[5]  # Related pair
  c[4].p1, c[4].p2 = p[6], p[7]
  c[5].p1, c[5].p2 = p[8], p[9]  # Unrelated pair
  c[6].p1, c[6].p2 = p[10], p[11]
  c[7].p1, c[7].p2 = p[10], p[11]  # Related pair

  mtr = gen.StockMater(incest_generations=1)
  mating_list, incest_list = mtr.mate_in_sequence(list(c))
  assert mating_list == [[c[5], c[4]], [c[1], c[0]]], mating_list
  assert incest_list == [[c[7], c[3]], [c[6], c[2]]], incest_list


def test_stock_breeder():
  """Stock breeding algorithm."""
  # These tests use random seeds. Need to figure out if there is a way to make these more deterministic
  brdr = gen.StockBreeder(child_factor=2.0, master_seed=3)
  p1 = gen.Sample(0, 0)
  p1.fitness = 0
  p2 = gen.Sample(0, 1)
  p2.fitness = 0

  c = brdr.breed([[p1, p2]], generation=1)
  assert len(c) == 2  # Perfectly average parents
  assert [c[0].p1, c[0].p2] == [p1, p2]

  p2.fitness = 1
  brdr = gen.StockBreeder(child_factor=2.0, master_seed=3)
  c = brdr.breed([[p1, p2]], generation=1)
  assert len(c) == 3

  p3 = gen.Sample(0, 2)
  p3.fitness = 0
  p4 = gen.Sample(0, 3)
  p4.fitness = 0

  c = brdr.breed([[p1, p2], [p3, p4]], generation=1)
  assert len(c) == 5
  assert c[0].generation == 1
  assert c[4].serial == 4
  assert [c[4].p1, c[4].p2] != [p1, p2]
  assert [c[4].p1, c[4].p2] == [p3, p4]


def test_stock_culler():
  """Stock culling algorithm."""
  p = [gen.Sample(0, n) for n in range(10)]
  p[3].fitness = -0.25
  p[9].fitness = -0.5
  p[5].fitness = -1.0

  cull = gen.StockCuller(8)
  cull.cull(p)

  assert len(p) == 8
  assert p[0].fitness == 0
  assert p[3].serial == 3
  assert p[5].serial == 6