import React from "react";

export default ({ handleChange, handleSubmit, post }) => {
  return (
    <div className="container">
      <div className="jumbotron jumbotron-fluid">
        <div
          className="form-container"
          style={{
            display: "flex",
            justifyContent: "center",
            alignItems: "center"
          }}
        >
          <form>
            <div className="form-group">
              <input
                className="col-sm-12"
                name="name"
                onChange={handleChange}
                type="text"
                value={post.name}
                placeholder="Gene Name or NCBI gene Id"
                width="75px"
              />
            </div>

            <button className="btn btn-primary" onClick={handleSubmit}>
              Submit
            </button>
          </form>
        </div>
      </div>
    </div>
  );
};
