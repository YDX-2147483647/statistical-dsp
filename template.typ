#import "@preview/oxifmt:0.2.0": strfmt

#let project(
  title: none,
  author: (name: none),
  date: none,
  //! Strings are passed to `strfmt`, and contents or other types are shown literally.
  info: (),
  body,
) = {
  // Document's basic properties
  set document(author: author.name, title: title)
  set page(numbering: "1 / 1", number-align: center)

  // Font families
  let body-font = ("Linux Libertine", "Source Han Serif")
  let sans-font = ("Inria Sans", "Source Han Sans")

  // Text formats
  set text(font: body-font, lang: "en")
  show heading: set text(font: sans-font)
  set heading(numbering: "1.1")

  //! Title page

  //! Headings
  align(left, image("assets/logo.png", height: 5em))
  // The logo was modified from images in BIThesis.
  // https://github.com/BITNP/BIThesis/tree/a4daadbb00ea8cee8587ba79f8ed79ddbc3da8a9/templates/undergraduate-thesis/images

  v(1fr, weak: true)
  align(
    center,
    text(font: sans-font, 2em, weight: "bold")[Undergraduate Experiment Report],
  )

  v(1fr, weak: true)
  align(center)[
    #set text(font: body-font, 2em, weight: "bold")

    Laboratory:

    #set text(style: "italic")

    #title
  ]

  //! Information table
  let formatted_info = info.pairs().map(((key, value))=>{
    if type(value) == str {
      (key, strfmt(value, date: date, ..author))
    } else {
      (key, value)
    }
  })

  let n_cols = 2
  let n_rows = calc.ceil(formatted_info.len() / n_cols)

  v(1fr, weak: true)
  align(center, table(
    columns: (auto, auto) * n_cols,
    align: center + horizon,
    inset: (x: 1em, y: 1em),
    ..range(n_rows).map(r=>{
      range(n_cols).map(c=>{
        let i = r + n_rows * c
        if i < formatted_info.len() {
          let pair = formatted_info.at(i)
          ({
            set text(weight: "bold")
            pair.at(0)
          }, pair.at(1))
        }
      })
    }).flatten(),
  ))

  v(2fr, weak: true)
  pagebreak()

  //! Main body

  set par(justify: true)

  outline(indent: 2em)

  body
}
