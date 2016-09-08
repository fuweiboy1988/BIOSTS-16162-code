#
# Figure 1 of supplemental material
#
## package
library("partykit")

## node/tree structure (without data)
## info has mathematical expression to be plotted
pn <- partynode(1L, split = partysplit(1L, breaks = 2), kids = list(
  partynode(2L, split = partysplit(2L, breaks = 1), kids = list(
    partynode(3L, info = expression(tilde(T)[1])),
    partynode(4L, info = expression(tilde(T)[2])))),
  partynode(5L, split = partysplit(3L, breaks = 1), kids = list(
    partynode(6L, info = expression(tilde(T)[3])),
    partynode(7L, info = expression(tilde(T)[4]))))))

## create party object: node/tree + data (empty here)
py <- party(node = pn,
  data = data.frame(X1 = numeric(0), X2 = numeric(0), X3 = numeric(0)))

pdf("c:/papers/trees/survival/biostatisticsrevision/tree_structure.pdf")

## tweaks: suppress IDs and show mathematical expressions directly
plot(py, ip_args = list(id = FALSE),
  tp_args = list(id = FALSE, FUN = identity, width = 7))

dev.off()

#
# Figure 21 of supplemental material
#


pdf("c:/papers/trees/survival/biostatisticsrevision/TV_tree_structure1.pdf")


pn1 <- partynode(1L, split = partysplit(1L, breaks = 0), kids = list(
  partynode(2L, split = partysplit(2L, breaks = 0), kids = list(
    partynode(3L, info = expression(tilde(T)[4])),
    partynode(4L, info = expression(tilde(T)[3])))),
  partynode(5L, split = partysplit(2L, breaks = 0), kids = list(
    partynode(6L, info = expression(tilde(T)[2])),
    partynode(7L, info = expression(tilde(T)[1]))))))

## create party object: node/tree + data (empty here)
py1 <- party(node = pn1,
  data = data.frame(X1 = numeric(0), X2 = numeric(0), X3 = numeric(0)))

plot(py1, ip_args = list(id = FALSE),
  tp_args = list(id = FALSE, FUN = identity, width = 7))

dev.off()

pdf("c:/papers/trees/survival/biostatisticsrevision/TV_tree_structure2.pdf")

pn2 <- partynode(1L, split = partysplit(2L, breaks = 0), kids = list(
  partynode(2L, split = partysplit(1L, breaks = 0), kids = list(
    partynode(3L, info = expression(tilde(T)[4])),
    partynode(4L, info = expression(tilde(T)[2])))),
  partynode(5L, split = partysplit(1L, breaks = 0), kids = list(
    partynode(6L, info = expression(tilde(T)[3])),
    partynode(7L, info = expression(tilde(T)[1]))))))

## create party object: node/tree + data (empty here)
py2 <- party(node = pn2,
  data = data.frame(X1 = numeric(0), X2 = numeric(0), X3 = numeric(0)))



plot(py2, ip_args = list(id = FALSE),
  tp_args = list(id = FALSE, FUN = identity, width = 7))

dev.off()

#
# Figure 22 of supplemental material
#


pdf("c:/papers/trees/survival/biostatisticsrevision/TV_tree_structure_ctn1.pdf")


pn3 <- partynode(1L, split = partysplit(1L, breaks = 0), kids = list(
  partynode(2L, split = partysplit(2L, breaks = 5), kids = list(
    partynode(3L, info = expression(tilde(T)[4])),
    partynode(4L, info = expression(tilde(T)[3])))),
  partynode(5L, split = partysplit(2L, breaks = 5), kids = list(
    partynode(6L, info = expression(tilde(T)[2])),
    partynode(7L, info = expression(tilde(T)[1]))))))

## create party object: node/tree + data (empty here)
py3 <- party(node = pn3,
  data = data.frame(X1 = numeric(0), X2 = numeric(0), X3 = numeric(0)))

plot(py3, ip_args = list(id = FALSE),
  tp_args = list(id = FALSE, FUN = identity, width = 7))

dev.off()

pdf("c:/papers/trees/survival/biostatisticsrevision/TV_tree_structure_ctn2.pdf")

pn4 <- partynode(1L, split = partysplit(2L, breaks = 5), kids = list(
  partynode(2L, split = partysplit(1L, breaks = 0), kids = list(
    partynode(3L, info = expression(tilde(T)[4])),
    partynode(4L, info = expression(tilde(T)[2])))),
  partynode(5L, split = partysplit(1L, breaks = 0), kids = list(
    partynode(6L, info = expression(tilde(T)[3])),
    partynode(7L, info = expression(tilde(T)[1]))))))

## create party object: node/tree + data (empty here)
py4 <- party(node = pn4,
  data = data.frame(X1 = numeric(0), X2 = numeric(0), X3 = numeric(0)))



plot(py4, ip_args = list(id = FALSE),
  tp_args = list(id = FALSE, FUN = identity, width = 7))

dev.off()


