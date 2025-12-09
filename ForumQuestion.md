# How to prevent GraphicsRow from compressing individual graphs?

I have  graphs created with `Graph` that display OK when shown individually, 
but when I put them in `GraphicsRow`, they get compressed and the vertex 
labels become illegible.

## Minimal Example

```mathematica
(* Create vertices  *)
vertices = {{}, {2}, {3}, {4}, {2, 3}, {2, 4}, {3, 4}, {2, 3, 4}};

(* Example edges *)
edges = {{2, 3} -> {2}, {2, 4} -> {4}};

(* Create labels with subscripts *)
labels = {
  {} -> "",
  {2} -> Placed[Subscript[i, 12], Above],
  {3} -> Placed[Subscript[i, 1], Above],
  {4} -> Placed[Subscript[i, 2], Above],
  {2, 3} -> Placed[Row[{Subscript[i, 12], Subscript[i, 1]}], Below],
  {2, 4} -> Placed[Row[{Subscript[i, 12], Subscript[i, 2]}], Below],
  {3, 4} -> Placed[Row[{Subscript[i, 1], Subscript[i, 2]}], Below],
  {2, 3, 4} -> Placed[Row[{Subscript[i, 12], Subscript[i, 1], Subscript[i, 2]}], Below]
};

(* Create graph with fixed size *)
g1 = Graph[vertices, edges,
  VertexLabels -> labels,
  VertexSize -> 0.4,
  GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Top},
  ImageSize -> {300, 300},
  ImagePadding -> {{40, 40}, {30, 30}}
];

(* This displays fine - labels are readable *)
g1

(* Create three copies for demonstration *)
g2 = g1;
g3 = g1;

(* THIS IS THE PROBLEM: graphs get compressed, labels overlap *)
GraphicsRow[{g1, g2, g3}, Spacings -> 20]
```

## The Problem

When I display a single graph, it appears at 300×300 pixels with readable labels. But `GraphicsRow[{g1, g2, g3}, Spacings -> 20]` compresses each graph significantly, making the vertex labels overlap and become illegible.

I tried:
- Removing `ImageSize -> Full` from GraphicsRow
- Increasing `Spacings`
- Setting explicit `ImageSize -> 900` in GraphicsRow
- Increasing `ImagePadding` in the individual graphs

**Question:** 
How can I make `GraphicsRow` display each graph  
without compression, so the labels remain readable?
