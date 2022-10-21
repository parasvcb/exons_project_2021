const dataIPR = [
  {
    accession: "IPR00001",
    name: "Thiolase-like",
    source_database: "interpro",
    type: "homologous_superfamily",
    locations: [{ fragments: [{ start: -50, end: 1 }] }],
    color: "#cacaca",
  },
  {
    accession: "IPR00001",
    name: "Thiolase-like",
    source_database: "interpro",
    type: "homologous_superfamily",
    locations: [{ fragments: [{ start: 1, end: 45 }] }],
    color: "#7acaca",
  },
  {
    accession: "SF00001",
    name: "Thiolase-like",
    source_database: "interpro",
    type: "homologous_superfamily",
    locations: [
      {
        fragments: [
          { shape: 'discontinuosEnd', start: 127, end: 264 },
          { shape: 'discontinuosStart', start: 279, end: 385 },
          { shape: 'discontinuos', start: 400, end: 466 },
        ],
      },
    ],
    color: "#7acaca",
  },
];
const signatures = [
  {
    accession: "SF00001",
    name: "Thiolase-like",
    source_database: "interpro",
    type: "homologous_superfamily",
    locations: [
      {
        fragments: [
          { shape: 'discontinuosEnd', start: 127, end: 264 },
          { shape: 'discontinuosStart', start: 279, end: 385 },
          { shape: 'discontinuos', start: 400, end: 466 },
        ],
      },
    ],
    color: "#7acaca",
  },
  {
    accession: "PF00001",
    name: "Thiolase-like",
    source_database: "interpro",
    type: "homologous_superfamily",
    locations: [
      { fragments: [{ start: 146, end: 266 }] },
      { fragments: [{ start: 346, end: 466 }] },
    ],
    color: "#ca7aca",
  },
  {
    accession: "CD0001",
    name: "Thiolase-like",
    source_database: "interpro",
    type: "homologous_superfamily",
    locations: [{ fragments: [{ start: 146, end: 456 }] }],
    color: "#caca7a",
    //tooltip: "CD0001",
    //values:[{"position":10,"value":20}],
  },
];
export { dataIPR, signatures};
