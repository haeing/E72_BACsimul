#ifndef AeroHit_hh
#define AeroHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"

class G4AttDef;
class G4AttValue;

class AeroHit : public G4VHit
{
public:
  AeroHit();
  AeroHit(G4int i, G4double t);
  AeroHit(const AeroHit &right);
  virtual ~AeroHit();

  const AeroHit& operator=(const AeroHit &right);
  G4bool operator==(const AeroHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;

  void SetEdep(G4double de) {fEdep = de;}
  void AddEdep(G4double de) {fEdep += de;}
  G4double GetEdep() const {return fEdep;}

  G4int GetID() const {return fId;}

  void SetTime(G4double dt) {fTime = dt;}
  G4double GetTime() const {return fTime;}

  //void IncPhotonCount() {fPhotons++;}
  //G4double GetPhotonCount() {return fPhotons;}
  //void ClearPhotonCount() {fPhotons=0.0;}

  //void SetPos(G4ThreeVector xyz) {fPos = xyz;}
  //G4ThreeVector GetPos() const {return fPos;}
  
private:
  G4int fId;
  G4double fEdep;
  G4double fTime;
  //G4double fPhotons;
  G4ThreeVector fPos;
};

using AeroHitsCollection = G4THitsCollection<AeroHit>;

extern G4ThreadLocal G4Allocator<AeroHit>* AeroHitAllocator;

inline void* AeroHit::operator new(size_t)
{
  if (!AeroHitAllocator){
    AeroHitAllocator = new G4Allocator<AeroHit>;
  }
  return (void*)AeroHitAllocator-> MallocSingle();
}

inline void AeroHit::operator delete(void* aHit)
{
  AeroHitAllocator->FreeSingle((AeroHit*) aHit);
}

#endif
  
