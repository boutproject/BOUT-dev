virtual const Field3D indexDDX(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method) override;
virtual const Field3D indexDDY(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method) override;
virtual const Field3D indexDDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               bool ignored) override;
virtual const Field2D indexDDX(const Field2D &f) override;
virtual const Field2D indexDDY(const Field2D &f) override;
virtual const Field3D indexD2DX2(const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) override;
virtual const Field3D indexD2DY2(const Field3D &f, CELL_LOC outloc,
                                 DIFF_METHOD method) override;
virtual const Field3D indexD2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                 bool ignored) override;
virtual const Field2D indexD2DX2(const Field2D &f) override;
virtual const Field2D indexD2DY2(const Field2D &f) override;
virtual const Field3D indexVDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field3D indexVDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field3D indexVDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                DIFF_METHOD method);
virtual const Field3D indexVDDZ(const Field &v, const Field &f, CELL_LOC outloc,
                                DIFF_METHOD method) override {
  return indexVDDZ(dynamic_cast<const Field3D &>(v), dynamic_cast<const Field3D &>(f),
                   outloc, method);
}
virtual const Field2D indexVDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field2D indexVDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field3D indexFDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field3D indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field3D indexFDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                DIFF_METHOD method) override;
virtual const Field2D indexFDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
virtual const Field2D indexFDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                                DIFF_METHOD method, REGION ignored) override;
