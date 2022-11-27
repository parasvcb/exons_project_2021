const dataIPR = [
  {
    accession: "IPR00001",
    locations: [{ fragments: [{ start: 1, end: 45 }] }],
    color: "#7acaca",
  },
  {
    accession: "IPR00002",
    locations: [{ fragments: [{ start: 112, end: 200 }] }],
    color: "#7acaca",
    opacity: 0.5
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
    opacity: "0.5"
  },
  {
    accession: "CD0001",
    locations: [{ fragments: [{ start: 100, end: 200 }] }],
    color: "yellow",
  },
];
const signatures = [

  {
    accession: "CD0001",
    locations: [{ fragments: [{ start: 146, end: 456 }] }],
    color: "#caca7a",
  },
  {
    accession: "CD0001",
    locations: [{ fragments: [{ start: 500, end: 700 }] }],
    color: "#caca7a",
  },
  {
    accession: "CD0001",
    locations: [{ fragments: [{ start: 130, end: 150 }] }],
    color: "#caca7a",
  },
  {
    accession: "CD0001",
    locations: [{ fragments: [{ start: 100, end: 120 }] }],
    color: "#caca7a",
  },
  {
    accession: "CD0001",
    locations: [{ fragments: [{ start: 100, end: 200 }] }],
    color: "yellow",
  },
];
export { dataIPR, signatures};
