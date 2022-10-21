import React, { Component } from "react";
import SeqExonMatching from "./UserTranscriptMathematics/SeqExonMatching";
import {
  Form,
  FormControl,
  FormGroup,
  ControlLabel,
  HelpBlock,
  Checkbox,
  Radio,
  Button
} from "react-bootstrap";

export class UserTranscripts extends Component {
  constructor(props) {
    super(props);
    this.state = {
      value: "",
      proceed: false,
      formseq: false
    };
    this.handleChange = this.handleChange.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
    this.seqFormat = this.seqFormat.bind(this);
  }

  handleChange(event) {
    this.setState({ value: event.target.value });
  }

  handleSubmit(event) {
    let sequformatted = "";
    console.log("event", event.value);
    event.preventDefault();
    this.setState({ proceed: true });
    sequformatted = this.seqFormat(this.state.value);
    this.setState({ formseq: sequformatted });
  }
  seqFormat(fasta) {
    let seq = "";
    let repat = /^>/g;
    fasta
      .split("\n")
      .map(ex => (seq = !ex.match(repat) ? seq + ex.trim() : ""));
    console.log(seq, seq.length);
    return seq;
  }

  render() {
    return this.state.proceed === false ? (
      <form onSubmit={this.handleSubmit}>
        <Form.Group controlId="exampleForm.ControlTextarea1">
          <Form.Label>Enter your sequence in FASTA format below</Form.Label>
          <Form.Control as="textarea" rows="10" onChange={this.handleChange} />
          <Button
            className="btn btn-primary btn-large centerButton"
            type="submit"
          >
            Send
          </Button>
        </Form.Group>
      </form>
    ) : (
      <>
        <SeqExonMatching
          fasta={this.state.formseq}
          transcripts={this.props.transet}
          exOn={this.props.exOn}
          exoncodes={this.props.exonCodes}
          exonAlerts={this.props.exonAlerts}
        />
      </>
    );
  }
}

export default UserTranscripts;
