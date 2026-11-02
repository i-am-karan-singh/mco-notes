#import "@preview/theorion:0.4.1": *
#import cosmos.clouds: *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#let lecturetitle(title, number) = {
  align(center)[
    #box(
      inset: 1em,
      stroke: black
    )[
      *47-851: Modern Convex Optimization* #h(1fr) Lecture #number

      #text(size: 1.5em)[
        #title
      ]

      #h(1fr) Lecturer: Karan Singh
    ]
    #v(1em)
  ]
}

#let lecture(title, number, doc) = {
  set page(paper: "us-letter")
  set heading(numbering: "1.")
  set text(font: "New Computer Modern")
  set par(justify: true)
  set outline(target:
  	selector(heading)
   	.after(label("start" + str(number)))
   	.before(label("end" + str(number)))
  )

  set-inherited-levels(0)

  show heading: set block(above: 1.4em, below: 1em)
  show link: underline
  show: show-theorion

  [#metadata(none)#label("start" + str(number))]
  lecturetitle(title, number)

  doc
  [#metadata(none)#label("end" + str(number))]
}

#let citet(key) = cite(key, form: "prose")
#let claim = proposition.with(supplement: [Claim])
#let angle(a, b) = $⟨#a, #b⟩$
#let BQP = $upright("BQP")$
#let SDP = $upright("SDP")$
#let OPT = $upright("OPT")$
