#if !defined  ROBUSTLMM_INDEX_H__
#define  ROBUSTLMM_INDEX_H__

#include "misc.h"

class BaseIndex;
class BlockTypeIndex;
class BlockIndex;
class RandomEffectIndex;
class ThetaIndex;
class IndexMapper;

class NonAssignable {
public:
  NonAssignable(NonAssignable const&) = delete;
  NonAssignable& operator=(NonAssignable const&) = delete;
  NonAssignable() {}
};

class BaseIndex : public NonAssignable {
  friend class RandomEffectIndex;
  friend class BlockIndex;
  friend class ThetaIndex;
  friend class BlockTypeIndex;
  friend class BaseIndexWithBlockType;

private:
  const unsigned index_;
  const IndexMapper* indexMapper_;

protected:
  BaseIndex(const unsigned index, const IndexMapper* indexMapper);
  const IndexMapper* getIndexMapper() const;

public:
  unsigned getIndex() const;
};

class BaseIndexWithBlockType : public BaseIndex {
  friend class RandomEffectIndex;
  friend class BlockIndex;
  friend class ThetaIndex;

private:
  const unsigned blockTypeIndex_;

protected:
  BaseIndexWithBlockType(const unsigned index, const unsigned blockTypeIndex,
                         const IndexMapper* indexMapper);

public:
  const BlockTypeIndex* getBlockTypeIndex() const;
  bool isDiagonalBlockType() const;
  bool isNonDiagonalBlockType() const;
  unsigned getBlockTypeDimension() const;
  bool isBlockTypeDropped() const;
};

class RandomEffectIndex : public BaseIndexWithBlockType {
  friend class IndexMapper;

private:
  const unsigned blockIndex_;

protected:
  RandomEffectIndex(const unsigned randomEffectIndex, const unsigned blockTypeIndex,
                    const IndexMapper* indexMapper, const unsigned blockIndex);

public:
  const BlockIndex* getBlockIndex() const;
};

class BlockIndex : public BaseIndexWithBlockType {
  friend class IndexMapper;

private:
  const std::vector<RandomEffectIndex*> randomEffects_;

protected:
  BlockIndex(const unsigned blockIndex, const unsigned blockTypeIndex,
             const IndexMapper* indexMapper,
             const std::vector<RandomEffectIndex*> randomEffects);

public:
  const std::vector<RandomEffectIndex*>& getRandomEffects() const;
};

class ThetaIndex : public BaseIndexWithBlockType {
  friend class IndexMapper;

protected:
  ThetaIndex(const unsigned thetaIndex, const unsigned blockTypeIndex,
             const IndexMapper* indexMapper);
};

class BlockTypeIndex : public BaseIndex {
  friend class IndexMapper;

public:
  const unsigned dim_;

private:
  bool blockTypeDropped_;
  const std::vector<RandomEffectIndex*> randomEffects_;
  const std::vector<BlockIndex*> blocks_;
  const std::vector<ThetaIndex*> thetas_;

protected:
  BlockTypeIndex(const unsigned blockTypeIndex, const IndexMapper* indexMapper,
                 const unsigned dim,
                 const std::vector<RandomEffectIndex*> randomEffects,
                 const std::vector<BlockIndex*> blocks,
                 const std::vector<ThetaIndex*> thetas);

public:
  bool isNonDiagonal() const;
  bool isDiagonal() const;
  bool isDropped() const;
  const std::vector<RandomEffectIndex*>& getRandomEffects() const;
  const std::vector<BlockIndex*>& getBlocks() const;
  const std::vector<ThetaIndex*>& getThetas() const;
  unsigned getNumberOfBlocks() const;

private:
  bool isDropped(const VectorXd& theta) const;

protected:
  void updateDropped(const VectorXd& theta);
};

class IndexMapper : public NonAssignable {
  friend class RandomEffectIndex;
  friend class BlockIndex;
  friend class BaseIndexWithBlockType;
  friend class ThetaIndex;
  friend class BlockTypeIndex;

private:
  const MiVec ind_, k_;
  const VectorXi bBlockMap_;

  bool anyBlockTypeNonDiagonal_;
  bool anyBlockTypeDropped_;
  bool allBlockTypesDropped_;

  std::vector<std::unique_ptr<RandomEffectIndex>> randomEffects_;
  std::vector<std::unique_ptr<BlockIndex>> blocks_;
  std::vector<std::unique_ptr<ThetaIndex>> thetas_;
  std::vector<std::unique_ptr<BlockTypeIndex>> blockTypes_;

public:
  IndexMapper(const Rcpp::List& args);

private:
  void initialiseRandomEffects();
  void initialiseBlocks(const std::vector<VectorXi>& blockBMap);
  void initialiseThetas(const std::vector<VectorXi>& blockThetaMap);
  void initialiseBlockTypes(const MiVec& dim, const std::vector<VectorXi>& idx,
                            const std::vector<VectorXi>& blockThetaMap);

public:
  const MiVec& getInd() const;
  const MiVec& getK() const;
  const VectorXi& getBBlockMap() const;

  unsigned getNumberOfRandomEffects() const;
  unsigned getNumberOfBlocks() const;
  unsigned getNumberOfThetas() const;
  unsigned getNumberOfBlockTypes() const;

  const std::vector<std::unique_ptr<RandomEffectIndex> >& getRandomEffects() const;
  const std::vector<std::unique_ptr<BlockIndex> >& getBlocks() const;
  const std::vector<std::unique_ptr<ThetaIndex> >& getThetas() const;
  const std::vector<std::unique_ptr<BlockTypeIndex> >& getBlockTypes() const;
  std::vector<std::unique_ptr<BlockTypeIndex> >& getBlockTypes();

  unsigned getMaxBlockTypeDimension() const;
  bool isAnyBlockTypeNonDiagonal() const;
  bool isAnyBlockTypeDropped() const;
  bool areAllBlockTypesDropped() const;

  void updateBlockTypesDropped(const VectorXd& theta);

protected:
  const RandomEffectIndex* getRandomEffect(unsigned randomEffectIndex) const;
  const BlockIndex* getBlockIndex(unsigned blockIndex) const;
  const ThetaIndex* getThetaIndex(unsigned thetaIndex) const;
  const BlockTypeIndex* getBlockTypeIndex(unsigned blockTypeIndex) const;
};

template<class T1>
bool isContiguousVector(const std::vector<T1*>& vector) {
  if (vector.size() <= 1) return true;
  for (unsigned i = 0, end = vector.size() - 1; i < end; i++) {
    if (vector[i]->getIndex() + 1 != vector[i+1]->getIndex()) return false;
  }
  return true;
}

#endif
