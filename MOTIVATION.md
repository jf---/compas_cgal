You're right that the *core idea* is close to evident once you have the medial axis — MAT gives you the pocket's skeleton, the clearance disk at every point, and a natural variable-stepover structure. Elber–Cohen–Drake basically said exactly that in 2005 (MATHSM), and Held's Voronoi machinery makes the geometry almost fall out. So why is it still a live research topic and a product battleground?

A few reasons the "little and evident" version doesn't survive contact:

1. **The engagement angle is the hard part, not the loop.** Drawing trochoidal loops along the MAT is easy. Guaranteeing the *tool-engagement angle* stays under a cap everywhere is not — it depends on the already-machined region, not just the local clearance. Held's 2025 contribution over MATHSM is precisely the exact θ computation and the inverse problem "place the next machining circle so max θ = θ_max," which he can only do by bisection with no closed form. That's the non-obvious core, and it's exactly the axis this project attacks with an exact certificate instead of ε-bisection.

2. **Corners, bottlenecks, and merges.** The evident version breaks at run-merges (two engaged arcs fusing → O(1) jump in engagement the smooth model misses — the exact hole caught in this very thread), at neck/bottleneck regions where θ_max must be *reduced* to avoid chip re-cutting/welding in hard metal, and at concave transitions where MRR spikes. The "if you know MAT" intuition doesn't cover the discontinuities.

3. **It's an economics topic, not just a geometry one.** iMachining, VoluMill, Fusion Adaptive are commercial because holding constant engagement lets you run 3–5× the feed and full flute depth — huge cycle-time and tool-life money. When the payoff is that large, the last 10% of engagement smoothness, corner behavior, and robustness across arbitrary pockets is worth serious R&D even if the first 90% is a weekend project.

4. **MAT itself is the hidden monster.** "If you know MAT" is load-bearing. Exact medial axis of a region with circular-arc boundaries (Held's VRONI/ArcVRONI) is genuinely hard robustness work — degeneracies, near-tangencies, topology. Most of the difficulty people wave away lives there.

5. **Continuity, feasibility, and the physical machine.** G¹/G² smoothness for high feed, machinability checks, tool-size suggestion, linking moves, the ~50% non-cutting air problem — each is a real sub-topic.

So: the *skeleton→loops* insight is indeed evident with MAT in hand. What keeps it a topic is turning that into a *guaranteed* engagement bound that survives corners, merges, bottlenecks, and arbitrary arc-boundary pockets — plus the fact that the money makes the hard 10% worth fighting over. This project's whole angle is that the guarantee should be *provable*, which is exactly the part that isn't evident.

